module DN2

    using LinearAlgebra

    # Točke in uteži Gaussove kvadrature za [-1, 1]
    const quad_points = [
        0.038772417506050821933,
        0.116084070675255208483,
        0.192697580701371099716,
        0.268152185007253681141,
        0.341994090825758473007,
        0.413779204371605001525,
        0.483075801686178712909,
        0.549467125095128202076,
        0.612553889667980237953,
        0.671956684614179548379,
        0.727318255189927103281,
        0.778305651426519387695,
        0.824612230833311663196,
        0.865959503212259503821,
        0.902098806968874296728,
        0.932812808278676533361,
        0.957916819213791655805,
        0.977259949983774262663,
        0.990726238699457006453,
        0.998237709710559200350
    ]

    const quad_weights = [
        0.077505947978424811264,
        0.077039818164247965588,
        0.076110361900626242372,
        0.074723169057968264200,
        0.072886582395804059061,
        0.070611647391286779696,
        0.067912045815233903826,
        0.064804013456601038075,
        0.061306242492928939167,
        0.057439769099391551367,
        0.053227846983936824355,
        0.048695807635072232061,
        0.043870908185673271992,
        0.038782167974472017640,
        0.033460195282547847393,
        0.027937006980023401099,
        0.022245849194166957262,
        0.016421058381907888713,
        0.010498284531152813615,
        0.004521277098533191258
    ]

    """
    Poišče presečišče daljic p in q.
    """
    function segment_intersection(p1, p2, q1, q2; eps=1e-10)
        dp = p2 .- p1
        dq = q2 .- q1
        d = [-dp[2], dp[1]]
        denom = dot(dq, d)
        if abs(denom) < eps
            return nothing  # Vzporedne ali kolinearne
        end
        t = dot(q1 .- p1, d) / denom
        u = dot(q1 .- p1, [-dq[2], dq[1]]) / dot(dp, [-dq[2], dq[1]])
        if 0 ≤ t ≤ 1 && 0 ≤ u ≤ 1
            return t, u
        end
        return nothing
    end

    """
    evaluate_bezier(ctrls::Vector{Vector{Float64}}, t::Float64)

    Izračuna vrednosti Bezierjeve krivulje, dane s ctrls v točki t.
    """
    function evaluate_bezier(ctrls::Vector{Vector{Float64}}, t::Float64)
        tmp = deepcopy(ctrls)
        n = length(ctrls)
        for k in 1:n-1
            for i in 1:n-k
                tmp[i] = (1 - t) * tmp[i] .+ t * tmp[i+1]
            end
        end
        return tmp[1]
    end

    """
    Izračuna točke za vzorčenje krivulje.
    """
    function bezier_polyline(ctrls::Vector{Vector{Float64}}, samples::Int = 1000)
        ts = range(0.0, 1.0; length=samples)
        points = [evaluate_bezier(ctrls, t) for t in ts]
        return ts, points
    end    
    
    """
    Poišče približno presečišče daljic na krivulji.
    """
    function find_precise_intersection(ctrls::Vector{Vector{Float64}}, samples::Int=1000)
        ts, pts = bezier_polyline(ctrls, samples)
        for i in 1:samples-2
            p1, p2 = pts[i], pts[i+1]
            for j in i+2:samples-1  # Preskoči sosednje ali prekrivajoče
                q1, q2 = pts[j], pts[j+1]
                result = segment_intersection(p1, p2, q1, q2)
                if result !== nothing
                    (u1, u2) = result
                    t1 = (1 - u1)*ts[i] + u1*ts[i+1]
                    t2 = (1 - u2)*ts[j] + u2*ts[j+1]
                    return t1, t2
                end
            end
        end
        return nothing
    end

    """
    Iz začetnih t1 in t2 z bisekcijo izboljša oceno presečišča.
    """
    function refine_intersection(ctrls::Vector{Vector{Float64}}, t1::Float64, t2::Float64;
                                tol::Float64 = 1e-10, max_iter::Int = 100)

        function refine_single(fixed_t, lo, hi, fixed_is_t1)
            for _ in 1:max_iter
                mid = (lo + hi) / 2
                p_fixed = evaluate_bezier(ctrls, fixed_t)
                p_mid   = evaluate_bezier(ctrls, mid)
                dist_lo = norm(p_fixed .- evaluate_bezier(ctrls, lo))
                dist_hi = norm(p_fixed .- evaluate_bezier(ctrls, hi))
                dist_mid = norm(p_fixed .- p_mid)

                if dist_lo < dist_hi
                    hi = mid
                else
                    lo = mid
                end

                if abs(hi - lo) < tol
                    return (lo + hi) / 2
                end
            end
            return (lo + hi) / 2
        end

        # Izboljšaj t1 pri fiksnem t2
        t1 = refine_single(t2, t1 - 0.01, t1 + 0.01, true)
        # Izboljšaj t2 pri fiksnem t1
        t2 = refine_single(t1, t2 - 0.01, t2 + 0.01, false)

        return t1, t2
    end

    """
    Izračuna predznačeno ploščino po Greeovem teoremu.
    """
    function signed_area(ctrls::Vector{Vector{Float64}}; a=0.0, b=1.0)
        n = length(ctrls) - 1

        # Derivative control points
        dx = [n * (ctrls[i+1][1] - ctrls[i][1]) for i in 1:n]
        dy = [n * (ctrls[i+1][2] - ctrls[i][2]) for i in 1:n]

        # Bézier evaluation (vector and scalar)
        function B(ctrls::Vector{Float64}, t)
            m = length(ctrls) - 1
            sum(binomial(m, i) * (1 - t)^(m - i) * t^i * ctrls[i+1] for i in 0:m)
        end

        function B_vec(ctrls::Vector{Vector{Float64}}, t)
            m = length(ctrls) - 1
            sum(binomial(m, i) * (1 - t)^(m - i) * t^i .* ctrls[i+1] for i in 0:m)
        end

        # The integrand for Green's theorem
        integrand(t) = begin
            pt = B_vec(ctrls, t)
            x, y = pt[1], pt[2]
            ẋ = B(dx, t)
            ẏ = B(dy, t)
            x * ẏ - y * ẋ
        end

        # Gauss–Legendre quadrature over [a, b]
        s = 0.0
        for (w, t) in zip(quad_weights, quad_points)
            t = a + (b - a) * t
            s += w * integrand(t)
        end

        return 0.5 * (b - a) * s
    end

    """
    Krovna funkcija za izračun ploščine zanke, dane s kontrolnimi točkami ctrls.
    """
    function compute_area_bezier(ctrls::Vector{Vector{Float64}})
        intersections = find_precise_intersection(ctrls)
        if intersections === nothing
            return signed_area(ctrls)
        else
            t1, t2 = intersections
            t1, t2 = refine_intersection(ctrls, t1, t2)
            if t1 > t2
                t1, t2 = t2, t1
            end
            return signed_area(ctrls; a=t1, b=t2)
        end
    end

    export signed_area, compute_area_bezier, find_self_intersections, evaluate_bezier

end 
