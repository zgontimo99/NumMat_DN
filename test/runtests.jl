using DN2, Test

@testset "Examples" begin
    pts = [[0.0, 0.0] for _ in 1:8]
    x_pts = [p[1] for p in pts]
    y_pts = [p[2] for p in pts]
    xpoly = bezier_polynomial(x_pts)
    ypoly = bezier_polynomial(y_pts)
    area_poly = cross_term(xpoly, ypoly)
    area = 0.5 * poly_integrate(area_poly)
    @test abs(area) < 1e-10
    pts = [[0.0, 0.0], [1.0, 1.0], [2.0, 3.0], [1.0, 4.0], [0.0, 4.0], [-1.0, 3.0], [0.0, 1.0], [1.0, 0.0]]
    x_pts = [p[1] for p in pts]
    y_pts = [p[2] for p in pts]
    xpoly = bezier_polynomial(x_pts)
    ypoly = bezier_polynomial(y_pts)
    area_poly = cross_term(xpoly, ypoly)
    area1 = 0.5 * poly_integrate(area_poly)
    pts = reverse(pts)
    x_pts = [p[1] for p in pts]
    y_pts = [p[2] for p in pts]
    xpoly = bezier_polynomial(x_pts)
    ypoly = bezier_polynomial(y_pts)
    area_poly = cross_term(xpoly, ypoly)
    area2 = 0.5 * poly_integrate(area_poly)
    @test abs(area1 + area2) < 1e-10
end

