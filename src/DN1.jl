module DN1

    using LinearAlgebra, Graphs, Plots

    struct RedkaMatrika{T}
        V::Matrix{T}      # Values
        I::Matrix{Int}    # Column indices
    end

    function build_indexed_matrix(A::Matrix)
        n, m = size(A)
        max_nnz_per_row = maximum(count(!iszero, A[i, :] ) for i in 1:n)
    
        V = zeros(eltype(A), n, max_nnz_per_row)
        I = zeros(Int, n, max_nnz_per_row)
        for i in 1:n
            k = 1
            for j in 1:m
                if A[i, j] != 0
                    V[i, k] = A[i, j]
                    I[i, k] = j
                    k += 1
                end
            end
        end
    
        return RedkaMatrika(V, I)
    end

    function conj_grad(A::RedkaMatrika, b::AbstractVector; tol=1e-8, maxiter=length(b))
        n = length(b)
        x = copy(b)
        r = b - A * x
        p = copy(r)
        rsold = dot(r, r)
    
        for i in 1:maxiter
            Ap = A * p
            alpha = rsold / dot(p, Ap)
            x += alpha * p
            r -= alpha * Ap
            rsnew = dot(r, r)
    
            if sqrt(rsnew) < tol
                return x
            end
    
            p = r + (rsnew / rsold) * p
            rsold = rsnew
        end
    
        return x  # If convergence not reached
    end

    import Base: getindex, setindex!, firstindex, lastindex, *

    function getindex(A::RedkaMatrika, i::Int, j::Int)
        row_indices = A.I[i, :]
        row_values = A.V[i, :]
        for k in eachindex(row_indices)
            if row_indices[k] == j
                return row_values[k]
            end
        end
        return zero(eltype(A.V))  # Implicit zero for not-stored entries
    end

    function setindex!(A::RedkaMatrika, val, i::Int, j::Int)
        row_indices = A.I[i, :]
        for k in eachindex(row_indices)
            if row_indices[k] == j
                A.V[i, k] = val
                return A
            end
        end
        error("Cannot insert new element at position ($i, $j): no slot reserved")
    end

    firstindex(A::RedkaMatrika) = (1, 1)
    lastindex(A::RedkaMatrika) = (size(A.V, 1), maximum(A.I))  # upper bound on columns

    function *(A::RedkaMatrika{T}, x::AbstractVector{T}) where T
        n = size(A.V, 1)
        result = zeros(T, n)
        for i in 1:n
            for k in 1:size(A.V, 2)
                col = A.I[i, k]
                if col == 0
                    continue  # Skip unused slots
                end
                result[i] += A.V[i, k] * x[col]
            end
        end
        return result
    end

    export RedkaMatrika, build_indexed_matrix, conj_grad
end