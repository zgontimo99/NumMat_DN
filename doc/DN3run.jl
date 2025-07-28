using DN3, Plots, Printf

u0 = [1.0, 0.0]  # začetna x in dx/dt
dt = 1e-5
tspan = (0.0, 100.0)  # Zagotovimo stabilen limitni cikel

t, u = rk4(van_der_pol!, u0, tspan, dt)

x = u[1, :]
v = u[2, :]

crossings = find_upward_crossings(t, x, v)

if length(crossings) ≥ 2
    period = crossings[end] - crossings[end - 1]
    @printf("Ocenjena perioda limitnega cikla ≈ %.10f\n", period)
else
    println("Premajhno število ciklov!")
end