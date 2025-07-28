using DN3, Test, Statistics

@testset "RK4: Harmonic Oscillator" begin
    function harmonic!(u, t)
        x, v = u
        return [v, -x]
    end

    u0 = [1.0, 0.0]
    t, u = rk4(harmonic!, u0, (0.0, 2π), 0.01)
    x = u[1, end]

    @test isapprox(x, 1.0; atol=1e-2)  # Should return to initial x
end

@testset "Zero-Crossing Interpolation" begin
    t = [0.0, 0.01]
    x = [-0.005, 0.005]
    v = [0.5, 0.5]

    crossings = find_upward_crossings(t, x, v)
    @test length(crossings) == 1
    @test isapprox(crossings[1], 0.005; atol=1e-6)
end

@testset "Stability Over Long Integration" begin
    t, u = rk4(van_der_pol!, [1.0, 0.0], (0.0, 1000.0), 1e-2)
    @test all(isfinite, u)
end

@testset "RK4 Time Step Sensitivity" begin
    periods = Float64[]
    for dt in [1e-3, 5e-4, 1e-4]
        t, u = rk4(van_der_pol!, [1.0, 0.0], (0.0, 100.0), dt)
        x, v = u[1, :], u[2, :]
        c = find_upward_crossings(t, x, v)
        p = c[end] - c[end-1]
        push!(periods, p)
    end
    @test maximum(periods) - minimum(periods) < 1e-4
end

@testset "Van der Pol: Convergence Test" begin
    u0 = [2.0, 0.0]
    t, u = rk4(van_der_pol!, u0, (0.0, 100.0), 1e-3)
    x = u[1, :]

    # Ensure solution is oscillatory in final segment
    minx = minimum(x[end-20000:end])
    maxx = maximum(x[end-20000:end])
    @test maxx - minx > 0.5  # Should be oscillating
end

@testset "Van der Pol: Edge Case Initial Conditions" begin
    for u0 in ([1e-10, 0.0], [20, 0.0], [1.0, 1e2])
        t, u = rk4(van_der_pol!, u0, (0.0, 50.0), 1e-3)
        x, v = u[1, end], u[2, end]
        @test isfinite(x)
        @test isfinite(v)
    end
end

@testset "Van der Pol: Period Accuracy" begin
    dt = 1e-4
    u0 = [1.0, 0.0]
    t, u = rk4(van_der_pol!, u0, (0.0, 100.0), dt)
    x, v = u[1, :], u[2, :]
    crossings = find_upward_crossings(t, x, v)

    @test length(crossings) > 2

    period = crossings[end] - crossings[end-1]

    reference = 10.203523625640202  # Known value for μ=4 from Julia's ODE solver
    @test isapprox(period, reference; atol=1e-6)
end
