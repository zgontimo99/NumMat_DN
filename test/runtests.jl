using DN3, Test, Statistics

# dx/dt = v, dv/dt = -x (simple harmonic oscillator)
function harmonic!(u, t)
    x, v = u
    return [v, -x]
end

@testset "RK4 Harmonic Oscillator Test" begin
    u0 = [1.0, 0.0]
    t, u = rk4(harmonic!, u0, (0.0, 2π), 0.01)
    x = u[1, end]
    @test isapprox(x, 1.0; atol=1e-2)
end

@testset "Van der Pol Convergence" begin
    u0 = [2.0, 0.0]
    t, u = rk4(van_der_pol!, u0, (0.0, 100.0), 1e-3)
    x = u[1, :]
    @test std(x[end-1000:end]) > 0.08  # Ensure it’s not flat
end