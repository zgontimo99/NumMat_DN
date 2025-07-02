module DN3

    using Plots
    μ = 0.0123
    # Define equations of motion for CR3BP in rotating frame
    function cr3bp_derivatives(u, t)
        x, y, z, vx, vy, vz = u

        r1 = sqrt((x + μ)^2 + y^2 + z^2)
        r2 = sqrt((x - 1 + μ)^2 + y^2 + z^2)

        ax = -(1 - μ)*(x + μ)/r1^3 - μ*(x - 1 + μ)/r2^3
        ay = -(1 - μ)*y/r1^3 - μ*y/r2^3
        az = -(1 - μ)*z/r1^3 - μ*z/r2^3

        return [vx, vy, vz, ax, ay, az]
    end

    function rk4(f, u0, tspan, dt)
        t0, tf = tspan
        N = Int(floor((tf - t0) / dt))
        t = range(t0, step=dt, length=N+1)
        
        u = zeros(length(u0), N+1)
        u[:, 1] = u0

        for i in 1:N
            k1 = f(u[:, i], t[i])
            k2 = f(u[:, i] .+ 0.5*k1*dt, t[i] + 0.5*dt)
            k3 = f(u[:, i] .+ 0.5*k2*dt, t[i] + 0.5*dt)
            k4 = f(u[:, i] .+ k3*dt, t[i] + dt)

            u[:, i+1] = u[:, i] .+ (dt/6)*(k1 .+ 2*k2 .+ 2*k3 .+ k4)
        end
        return t, u
    end


    export cr3bp_derivatives, rk4

end 
