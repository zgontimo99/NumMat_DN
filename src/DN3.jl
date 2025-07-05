module DN3

    using Plots

    function van_der_pol!(u, t)
        x, v = u
        dxdt = v
        dvdt = 4 * (1 - x^2) * v - x
        return [dxdt, dvdt]
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


    export rk4, van_der_pol!

end 
