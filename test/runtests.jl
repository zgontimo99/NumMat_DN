using DN1, Test, LinearAlgebra

@testset "IndexedMatrix Construction and Access" begin
    A = [4.0 1.0 0.0; 1.0 3.0 0.0; 0.0 0.0 2.0]
    indexed = build_indexed_matrix(A)

    @test size(indexed.V, 1) == 3
    @test indexed[1, 1] == 4.0
    @test indexed[1, 2] == 1.0
    @test indexed[1, 3] == 0.0  # Not stored → should return 0
    @test indexed[3, 3] == 2.0
end

@testset "IndexedMatrix Modification" begin
    A = [2.0 0.0 3.0; 0.0 5.0 0.0; 1.0 0.0 4.0]
    indexed = build_indexed_matrix(A)

    @test indexed[1, 3] == 3.0
    indexed[1, 3] = 9.0
    @test indexed[1, 3] == 9.0

    @test_throws ErrorException indexed[1, 2] = 5.0  # index not stored
end

@testset "Matrix-Vector Multiplication" begin
    A = [4.0 1.0 0.0; 1.0 3.0 0.0; 0.0 0.0 2.0]
    indexed = build_indexed_matrix(A)
    x = [1.0, 2.0, 3.0]
    result = indexed * x
    expected = A * x
    @test isapprox(result, expected, atol=1e-8)
end

@testset "Conjugate Gradient Solver" begin
    A = [4.0 1.0 0.0; 1.0 3.0 0.0; 0.0 0.0 2.0]
    b = [1.0, 2.0, 3.0]
    indexed = build_indexed_matrix(A)
    x = conj_grad(indexed, b)
    @test isapprox(A * x, b, atol=1e-8)
end

@testset "Edge Cases" begin
    # All-zero matrix
    A = zeros(3, 3)
    indexed = build_indexed_matrix(A)
    x = [1.0, 1.0, 1.0]
    result = indexed * x
    @test result == zeros(3)

    # Identity matrix
    A = [1.0 0 0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    b = [3.0, 2.0, 1.0]
    indexed = build_indexed_matrix(A)
    x = conj_grad(indexed, b)
    @test isapprox(x, b, atol=1e-8)

    # Matrix with only one nonzero per row
    A = diagm([1.0, 2.0, 3.0])
    indexed = build_indexed_matrix(A)
    b = [1.0, 2.0, 3.0]
    x = conj_grad(indexed, b)
    @test isapprox(x, [1.0, 1.0, 1.0], atol=1e-8)
end
