module DN2

    using LinearAlgebra

    # Binomial coefficient
    function binom(n, k)
        factorial(n) ÷ (factorial(k) * factorial(n - k))
    end

    # Multiply two polynomials (represented as coefficient arrays)
    function poly_mul(p, q)
        result = zeros(Float64, length(p) + length(q) - 1)
        for i in 1:length(p)
            for j in 1:length(q)
                result[i + j - 1] += p[i] * q[j]
            end
        end
        return result
    end

    # Compute the polynomial representation of a Bézier component (x or y)
    function bezier_polynomial(control_points)
        n = length(control_points) - 1
        poly = zeros(Float64, n + 1)
        for i in 0:n
            coeff = binom(n, i)
            # Bernstein basis: coeff * (1 - t)^(n - i) * t^i
            basis = poly_mul([coeff], poly_pow([1.0, -1.0], n - i))  # (1 - t)^(n - i)
            basis = poly_mul(basis, poly_pow([0.0, 1.0], i))         # t^i
            basis .= basis .* control_points[i + 1]
            poly = poly_add(poly, basis)
        end
        return poly
    end

    # Add two polynomials
    function poly_add(p, q)
        len = max(length(p), length(q))
        r = zeros(Float64, len)
        for i in 1:length(p)
            r[i] += p[i]
        end
        for i in 1:length(q)
            r[i] += q[i]
        end
        return r
    end

    # Raise a polynomial to a power
    function poly_pow(p, n)
        result = [1.0]
        for _ in 1:n
            result = poly_mul(result, p)
        end
        return result
    end

    # Derivative of a polynomial
    function poly_deriv(p)
        n = length(p)
        return [p[i] * (i - 1) for i in 2:n]
    end

    # Multiply two Bézier component polynomials and subtract
    function cross_term(xpoly, ypoly)
        dx = poly_deriv(xpoly)
        dy = poly_deriv(ypoly)
        term1 = poly_mul(xpoly, dy)
        term2 = poly_mul(ypoly, dx)
        return poly_add(term1, -1 .* term2)
    end

    # Integrate polynomial from 0 to 1
    function poly_integrate(p)
        s = 0.0
        for i in 1:length(p)
            s += p[i] / i
        end
        return s
    end

    export binom, poly_mul, bezier_polynomial, poly_add, poly_pow, poly_deriv, cross_term, poly_integrate

    # Control points
    control_pts = [
        [0.0, 0.0],
        [1.0, 1.0],
        [2.0, 3.0],
        [1.0, 4.0],
        [0.0, 4.0],
        [-1.0, 3.0],
        [0.0, 1.0],
        [1.0, 0.0]
    ]

    # Separate x and y control points
    x_pts = [p[1] for p in control_pts]
    y_pts = [p[2] for p in control_pts]

    # Get Bézier polynomials for x(t) and y(t)
    xpoly = bezier_polynomial(x_pts)
    ypoly = bezier_polynomial(y_pts)

    # Compute integrand x(t)*y'(t) - y(t)*x'(t)
    area_poly = cross_term(xpoly, ypoly)

    # Integrate from 0 to 1
    area = 0.5 * poly_integrate(area_poly)
    print(area)

end 
