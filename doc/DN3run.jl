using DN3, Plots, Printf

u0 = [1.0, 0.0]  # initial x and dx/dt
dt = 1e-4
tspan = (0.0, 100.0)  # Long enough to reach steady limit cycle

t, u = rk4(van_der_pol!, u0, tspan, dt)

x = u[1, :]
v = u[2, :]

function find_upward_crossings(t, x, v)
    crossings = Float64[]

    for i in 2:lastindex(x)
        if x[i-1] < 0 && x[i] >= 0 && v[i] > 0
            # Linear interpolation for better estimate:
            t0, t1 = t[i-1], t[i]
            x0, x1 = x[i-1], x[i]
            # t_cross = t0 - x0 * (t1 - t0) / (x1 - x0)
            t_cross = t0 + (0 - x0) * (t1 - t0) / (x1 - x0)
            push!(crossings, t_cross)
        end
    end
    return crossings
end

crossings = find_upward_crossings(t, x, v)

if length(crossings) ≥ 2
    period = crossings[end] - crossings[end - 1]
    @printf("Estimated limit cycle period ≈ %.10f\n", period)
else
    println("Not enough crossings to estimate period.")
end