module DN2

    using LinearAlgebra

    """Izračun binomskega koeficienta za n, k."""
    function binom(n, k)
        factorial(n) ÷ (factorial(k) * factorial(n - k))
    end

    """Polinomsko množenje polinomov p in q, predstavljenih kot polje koeficientov."""
    function poly_mul(p, q)
        result = zeros(Float64, length(p) + length(q) - 1)
        for i in eachindex(p)
            for j in eachindex(q)
                result[i + j - 1] += p[i] * q[j]
            end
        end
        return result
    end

    """Izračuna polinomsko predstavitev Bezierjeve komponente (x ali y)."""
    function bezier_polynomial(control_points)
        n = length(control_points) - 1
        result = Float64[]  # prilagodljiva velikost polinoma

        for i in 0:n
            coeff = binom(n, i)
            # Bernsteinova baza: coeff * (1 - t)^(n - i) * t^i
            one_minus_t = [1.0, -1.0]    # (1 - t)
            t_poly = [0.0, 1.0]          # t

            basis = poly_pow(one_minus_t, n - i)
            basis = poly_mul(basis, poly_pow(t_poly, i))
            basis = coeff .* basis
            basis = control_points[i + 1] .* basis

            result = poly_add(result, basis)
        end

        return result
    end

    """Izvede polinomsko seštevanje polinomov p in q."""
    function poly_add(p, q)
        len = max(length(p), length(q))
        r = zeros(Float64, len)
        for i in eachindex(p)
            r[i] += p[i]
        end
        for i in eachindex(q)
            r[i] += q[i]
        end
        return r
    end

    """Potencira polinom p na stopnjo n."""
    function poly_pow(p, n)
        result = [1.0]
        for _ in 1:n
            result = poly_mul(result, p)
        end
        return result
    end

    """Izračuna odvod polinoma p."""
    function poly_deriv(p)
        n = length(p)
        return [p[i] * (i - 1) for i in 2:n]
    end

    """Zmnnoži Bezierjeva polinoma za x in y, in nato odšteje enega od drugega."""
    function cross_term(xpoly, ypoly)
        dx = poly_deriv(xpoly)
        dy = poly_deriv(ypoly)
        term1 = poly_mul(xpoly, dy)
        term2 = poly_mul(ypoly, dx)
        return poly_add(term1, -1 .* term2)
    end

    """Izračuna integral polinoma med 0 in 1."""
    function poly_integrate(p)
        s = 0.0
        for i in eachindex(p)
            s += p[i] / i
        end
        return s
    end

    function compute_area(control_pts)
        # Ločitev x in y kontrolnih točk
        x_pts = [p[1] for p in control_pts]
        y_pts = [p[2] for p in control_pts]

        xpoly = bezier_polynomial(x_pts)
        ypoly = bezier_polynomial(y_pts)

        # Izračunaj integrand x(t)*y'(t) - y(t)*x'(t)
        area_poly = cross_term(xpoly, ypoly)

        # Intergriraj in deli z 2
        area = 0.5 * poly_integrate(area_poly)
        return area
    end

    export binom, poly_mul, bezier_polynomial, poly_add, poly_pow, poly_deriv, cross_term, poly_integrate, compute_area

end 
