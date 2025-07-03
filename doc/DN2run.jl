using DN2

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