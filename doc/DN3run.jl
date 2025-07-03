using DN3, Plots

tspan = (0.0, 4.0)

#u0 = [-0.1123, 0.0, 0.0, 0.0, 4.2175, 0.0]
#u0 = [-0.05615, 0.0, 0.0, 0.0, 6.5503, 0.0]
u0 = [-0.024, 0.0, 0.0, 0.0, 12.98, 0.0]

dt = 0.001

μ = 0.0123

t, sol = rk4(cr3bp_derivatives, u0, tspan, dt)

x = sol[1, :]
y = sol[2, :]

x_earth = -μ
x_moon  = 1 - μ

# Set axis limits (e.g., a bit beyond the Earth-Moon distance)
xlims = (-1.5, 1.5)
ylims = (-1.5, 1.5)

# Create the plot
plot(x, y,
    xlabel="x", ylabel="y",
    title="CR3BP Trajectory (RK4)",
    legend=false,
    lw=2,
    aspect_ratio=1,
    xlims=xlims,
    ylims=ylims,
    size=(600,600)
)

# Add Earth and Moon as fixed markers
scatter!([x_earth, x_moon], [0, 0],
    color=[:blue :gray],
    markersize=1,
    label=["Earth" "Moon"]
)

