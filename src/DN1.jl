module DN1

    using LinearAlgebra, Graphs

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
        x = zeros(eltype(b), n)
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
                return x, i
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
        for k in 1:length(row_indices)
            if row_indices[k] == j
                return row_values[k]
            end
        end
        return zero(eltype(A.V))  # Implicit zero for not-stored entries
    end

    function setindex!(A::RedkaMatrika, val, i::Int, j::Int)
        row_indices = A.I[i, :]
        for k in 1:length(row_indices)
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

    """
    G = krožna_lestev(n)
    Ustvari graf krožna lestev z `2n` točkami.
    """
    function krožna_lestev(n)
        G = SimpleGraph(2 * n)
        # prvi cikel
        for i = 1:n-1
            add_edge!(G, i, i + 1)
        end
        add_edge!(G, 1, n)
        # drugi cikel
        for i = n+1:2n-1
            add_edge!(G, i, i + 1)
        end
        add_edge!(G, n + 1, 2n)
        # povezave med obema cikloma
        for i = 1:n
            add_edge!(G, i, i + n)
        end
        return G
    end


    """
    A = matrika(G::AbstractGraph, sprem)
    Poišči matriko sistema linearnih enačb za vložitev grafa `G` s fizikalno metodo.
    Argument `sprem` je vektor vozlišč grafa, ki nimajo določenih koordinat.
    Indeksi v matriki `A` ustrezajo vozliščem v istem vrstnem redu,
    kot nastopajo v argumentu `sprem`.
    """
    function matrika(G, sprem)
        # preslikava med vozlišči in indeksi v matriki
        v_to_i = Dict([sprem[i] => i for i in eachindex(sprem)])
        m = length(sprem)
        A = zeros(m, m)
        for i = 1:m
            vertex = sprem[i]
            sosedi = neighbors(G, vertex)
            for vertex2 in sosedi
                if haskey(v_to_i, vertex2)
                    j = v_to_i[vertex2]
                    A[i, j] = 1
                end
            end
            A[i, i] = -length(sosedi)
        end
        return A
    end

    """
    b = desne_strani(G::AbstractGraph, sprem, koordinate)
    Poišči desne strani sistema linearnih enačb za eno koordinato vložitve grafa `G`
    s fizikalno metodo. Argument `sprem` je vektor vozlišč grafa, ki nimajo
    določenih koordinat. Argument `koordinate` vsebuje eno koordinato za vsa
    vozlišča grafa. Metoda uporabi le koordinato vozlišč, ki so pritrjena.
    Indeksi v vektorju `b` ustrezajo vozliščem v istem vrstnem redu,
    kot nastopajo v argumentu `sprem`.
    """
    function desne_strani(G, sprem, koordinate)
        set = Set(sprem)
        m = length(sprem)
        b = zeros(m)
        for i = 1:m
            v = sprem[i]
            for v2 in neighbors(G, v)
                if !(v2 in set) # dodamo le točke, ki so fiksirane
                    b[i] -= koordinate[v2]
                end
            end
        end
        return b
    end

    """
    vloži!(G::AbstractGraph, fix, točke)
    Poišči vložitev grafa `G` v prostor s fizikalno metodo. Argument `fix` vsebuje
    vektor vozlišč grafa, ki imajo določene koordinate. Argument `točke` je
    začetna vložitev grafa. Koordinate vozlišč, ki niso pritrjena, bodo nadomeščene
    z novimi koordinatami.
    Metoda ne vrne ničesar, ampak zapiše izračunane koordinate v matriko `točke`.
    """
    function vloži!(G, fix, točke)
        sprem = setdiff(vertices(G), fix)
        dim, _ = size(točke)
        A = build_indexed_matrix(matrika(G, sprem))
        for k = 1:dim
            b = desne_strani(G, sprem, točke[k, :])
            x = conj_grad(A, b) # matrika A je negativno definitna
            točke[k, sprem] = x
        end
    end


    export RedkaMatrika, build_indexed_matrix, conj_grad, krožna_lestev, matrika, desne_strani, vloži!
end