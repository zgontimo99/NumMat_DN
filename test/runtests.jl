using DN2, Test

# --- Helper: Close to zero ---
isapprox0(x; atol=1e-10) = isapprox(x, 0.0; atol=atol)

@testset "Bézier Curve Area Tests" begin

    @testset "Zero Area (degenerate cases)" begin
        # All points the same
        pts = [[0.0, 0.0] for _ in 1:8]
        @test isapprox0(compute_area_bezier(pts))

        # Straight line (collinear)
        pts = [[0.0, 0.0], [1.0, 1.0], [2.0, 2.0], [3.0, 3.0],
               [4.0, 4.0], [5.0, 5.0], [6.0, 6.0], [7.0, 7.0]]
        @test isapprox0(compute_area_bezier(pts))
    end

    @testset "Reversal flips sign" begin
        pts = [[0.0, 0.0], [1.0, 1.0], [2.0, 3.0], [1.0, 4.0],
               [0.0, 4.0], [-1.0, 3.0], [0.0, 1.0], [1.0, 0.0]]
        area1 = compute_area_bezier(pts)
        area2 = compute_area_bezier(reverse(pts))
        @test area1 > 0
        @test area2 < 0
    end

    @testset "Symmetry preservation" begin
        base_pts = [[Float64(x), Float64(y)] for (x, y) in
                    [(0,0), (1,1), (2,0), (1,-1), (0,-2), (-1,-1), (-2,0), (-1,1)]]

        @assert eltype(base_pts) == Vector{Float64}

        area1 = compute_area_bezier(base_pts)

        # Reflect across Y-axis with Float64 type enforcement
        pts_flipped = [[-p[1], p[2]] for p in base_pts]
        @assert eltype(pts_flipped) == Vector{Float64}

        area2 = compute_area_bezier(pts_flipped) * -1

        @test isapprox(area1, area2; atol=1e-10)
    end

    @testset "Known shape approximation" begin
        # Half circle-ish shape
        pts = [[1.0, 0.0], [1.0, 1.0], [-1.0, 1.0], [-1.0, 0.0]]
        area = compute_area_bezier(pts)
        true_area = 1.2  # approximation of semicircle area  ~ π/2
        @test isapprox(area, true_area; rtol=1e-2)  # generous due to shape approx
    end

    @testset "Randomized closed loop consistency" begin
        using Random
        Random.seed!(123)
        for _ in 1:5
            pts = [randn(2) for _ in 1:7]
            push!(pts, pts[1])  # close the loop
            area = compute_area_bezier(pts)
            @test isfinite(area)
        end
    end

    @testset "Different curve degrees" begin
        # Degree 2 (quadratic)
        pts2 = [[0.0, 0.0], [1.0, 2.0], [0.0, 4.0]]
        push!(pts2, pts2[1])
        while length(pts2) < 8
            push!(pts2, pts2[end])  # pad with duplicates
        end
        area = compute_area_bezier(pts2)
        @test isfinite(area)
    end
end

