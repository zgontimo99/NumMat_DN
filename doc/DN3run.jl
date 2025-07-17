using DN3, Plots, Printf

u0 = [1.0, 0.0]  # initial x and dx/dt
dt = 1e-5
tspan = (0.0, 100.0)  # Long enough to reach steady limit cycle

t, u = rk4(van_der_pol!, u0, tspan, dt)

x = u[1, :]
v = u[2, :]

crossings = find_upward_crossings(t, x, v)

if length(crossings) ≥ 2
    period = crossings[end] - crossings[end - 1]
    @printf("Estimated limit cycle period ≈ %.10f\n", period)
else
    println("Not enough crossings to estimate period.")
end