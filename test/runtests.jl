using DN1, Test

@testset "getindex" begin
    B = [2. 0;0 1]
    A = build_indexed_matrix(B)
    c = getindex(A, 2, 1)
    @test c ≈ 0.0
end

@testset "setindex" begin
    B = [2. 0;0 1]
    A = build_indexed_matrix(B)
    B = [2. 0;0 3]
    setindex!(A, 3.0, 2, 2)
    @test A[2, 2] ≈ B[2, 2]
end

@testset "Konjugirani gradienti" begin
    @test conj_grad(build_indexed_matrix([10. 0 2; 0 5 0; 1 0 3]),[16.0, 10.0, 10.0]) ≈ [1.1433910986649312, 1.873101830547942, 3.366340104101086]
end
