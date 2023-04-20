#!  Array comprehensions are wrong!
#!  y should be inner, and when 
#!  looping over ys, they should be reversed
if occursin("dennishb", homedir())
    using Pkg
    Pkg.activate("ComputationalPhysics", shared=true)
end
using LinearAlgebra
using Symbolics
using GLMakie; Makie.inline!(true)
update_theme!(resolution=(800, 800))
##
𝒢_Φ(x, xᶥ, y, yᶥ) = -1/4π * log((x-xᶥ)^2 + (y-yᶥ)^2)

let 
    @variables x xᶥ y yᶥ
    Dx = Differential(x)
    Dy = Differential(y)
    𝒢_Φ_expr = 𝒢_Φ(x, xᶥ, y, yᶥ)
    @variables x xᶥ y yᶥ
    Dx = Differential(x)
    Dy = Differential(y)

    𝒢_Ex_expr = -Dx(𝒢_Φ_expr) |> expand_derivatives
    𝒢_Ex_temp = build_function(𝒢_Ex_expr, [x xᶥ y yᶥ])|>eval
    global 𝒢_Ex(x, xᶥ, y, yᶥ) = 𝒢_Ex_temp([x, xᶥ, y, yᶥ])

    𝒢_Ey_expr = -Dy(𝒢_Φ_expr) |> expand_derivatives
    𝒢_Ey_temp = build_function(𝒢_Ey_expr, [x xᶥ y yᶥ])|>eval
    global 𝒢_Ey(x, xᶥ, y, yᶥ) = 𝒢_Ey_temp([x, xᶥ, y, yᶥ])
    # blay = -Dy(greenfunc) |> expand_derivatives
    # blax, blay
end

function parameters(;N_mesh, N_points=10, a=1, b=1)
    xmin, xmax = -2, 2
    xs = range(xmin, xmax, N_mesh)
    ymin, ymax = -2, 2
    ys = range(ymin, ymax, N_mesh)
    points = [[cos(c), sin(c)] for c in range(0, 2π, N_points+1)[begin:end-1]]
    return (;N_mesh, N_points, xs, ys, points, a, b)
end

##¤ a1) Plot the potential, E-field intensity, and E-field components
let
    p = parameters(N_mesh=100)
    (;N_mesh, xs, ys) = p
    xᶥ = 1
    yᶥ = 0
    𝒢s = [𝒢_Φ(x, xᶥ, y, yᶥ) for x in xs, y in ys]
    𝒢_Exs = [𝒢_Ex(x, xᶥ, y, yᶥ) for x in xs, y in ys]
    𝒢_Eys = [𝒢_Ey(x, xᶥ, y, yᶥ) for x in xs, y in ys]
    𝒢_Es_mag = hypot.(𝒢_Exs, 𝒢_Eys)
    fig = Figure()
    ax, hm = heatmap(fig[1, 1], xs, ys, 𝒢s, axis=(title="Φ", aspect=1, xlabel="x", ylabel="y"))
    Colorbar(fig[1, 2], hm)
    ax2, hm2 = heatmap(fig[1, 3], xs, ys, 𝒢_Es_mag, axis=(title="|E|", aspect=1, xlabel="x", ylabel="y"))
    Colorbar(fig[1, 4], hm2)
    ax3, hm3 = heatmap(fig[2, 1], xs, ys, 𝒢_Exs, axis=(title="Ex", aspect=1, xlabel="x", ylabel="y"))
    Colorbar(fig[2, 2], hm3)
    ax4, hm4 = heatmap(fig[2, 3], xs, ys, 𝒢_Eys, axis=(title="Ey", aspect=1, xlabel="x", ylabel="y"))
    Colorbar(fig[2, 4], hm4)
    colsize!(fig.layout, 1, Aspect(1, 1))
    colsize!(fig.layout, 3, Aspect(1, 1))
    resize_to_layout!(fig)
    display(fig)
end

##¤ a2)

Φ_ext(x, y, p) = p.a*x + p.b*y
E_ext_x(x, p) = -p.a
E_ext_y(y, p) = -p.b

##¤ b)
function norm_to_one(v::AbstractVector)
    most_extreme_val = v[findmax(abs, v)[2]]
    v ./= most_extreme_val
end
let 
    p = parameters(N_mesh=100, b=0)
    (;N_mesh, N_points, xs, ys, points) = p
    point_xs = first.(points)
    point_ys = last.(points)
    b⃗ = [-Φ_ext(point_xs[i], point_ys[i], p) for i in eachindex(points)]
    M = [i==j ? 0.0 : 𝒢_Φ(point_xs[i], point_xs[j], point_ys[i], point_ys[j]) for i in eachindex(points), j in eachindex(points)]# .|> ComplexF64
    M = issymmetric(M) ? Symmetric(M) : M
    ρs = M \ b⃗
    Φ_tot(x, y, ρs) = Φ_ext(x, y, p) + sum(ρs[i]*𝒢_Φ(x, point_xs[i], y, point_ys[i]) for i in eachindex(points))
    E_tot_x(x, y, ρs) = E_ext_x(x, p) + sum(ρs[i]*𝒢_Ex(x, point_xs[i], y, point_ys[i]) for i in eachindex(points))
    E_tot_y(x, y, ρs) = E_ext_y(y, p) + sum(ρs[i]*𝒢_Ey(x, point_xs[i], y, point_ys[i]) for i in eachindex(points))
    E_tot_mag(x, y, ρs) = hypot(E_tot_x(x, y, ρs), E_tot_y(x, y, ρs))

    fig = Figure()
    ax = Axis(fig[1, 1], aspect=1)
    Φ_tots = [Φ_tot(x, y, ρs) for x in xs, y in ys]
    hm = heatmap!(ax, xs, ys, Φ_tots)
    plt = scatter!(ax, Point2f.(points), marker=[ρ≥0 ? '+' : '-' for ρ in ρs], markersize=40 .* abs.(norm_to_one(ρs)), color=[ρ≥0 ? :red : :blue for ρ in ρs])
    Colorbar(fig[1, 2], hm, label="Φ_tot")

    E_tot_mags = [E_tot_mag(x, y, ρs) for x in xs, y in ys]
    ax2 = Axis(fig[1, 3], aspect=1)
    hm2 = heatmap!(ax2, xs, ys, E_tot_mags)#, colorrange=(0, 3))
    Colorbar(fig[1, 4], hm2, label="|E_tot|")

    E_tot_xs = [E_tot_x(x, y, ρs) for x in xs, y in ys]
    ax3 = Axis(fig[2, 1], aspect=1)
    hm3 = heatmap!(ax3, xs, ys, E_tot_xs)#, colorrange=(-3, 3))
    Colorbar(fig[2, 2], hm3, label="E_tot_x")

    E_tot_ys = [E_tot_y(x, y, ρs) for x in xs, y in ys]
    ax4 = Axis(fig[2, 3], aspect=1)
    hm4 = heatmap!(ax4, xs, ys, E_tot_ys)#, colorrange=(-3, 3))
    Colorbar(fig[2, 4], hm4, label="E_tot_y")

    # plt = scatter!(Point2f.(points), marker=[ρ≥0 ? '+' : '-' for ρ in ρs], markersize=40 .* abs.(norm_to_one(ρs)), color=[ρ≥0 ? :red : :blue for ρ in ρs])
    
    colsize!(fig.layout, 1, Aspect(1, 1))
    colsize!(fig.layout, 3, Aspect(1, 1))
    resize_to_layout!(fig)
    display(fig)
end

##! Investigations:
let 
    p = parameters(N_mesh=100, b=0)
    (;N_mesh, N_points, xs, ys, points) = p
    point_xs = first.(points)
    point_ys = last.(points)
    b⃗ = [-Φ_ext(point_xs[i], point_ys[i], p) for i in eachindex(points)]
    M = [i==j ? 0.0 : 𝒢_Φ(point_xs[i], point_xs[j], point_ys[i], point_ys[j]) for i in eachindex(points), j in eachindex(points)]# .|> ComplexF64
    M = issymmetric(M) ? Symmetric(M) : M
    ρs = M \ b⃗
    Φ_tot(x, y, ρs) = Φ_ext(x, y, p) + sum(ρs[i]*𝒢_Φ(x, point_xs[i], y, point_ys[i]) for i in eachindex(points))
    E_tot_x(x, y, ρs) = E_ext_x(x, p) + sum(ρs[i]*𝒢_Ex(x, point_xs[i], y, point_ys[i]) for i in eachindex(points))
    E_tot_y(x, y, ρs) = E_ext_y(y, p) + sum(ρs[i]*𝒢_Ey(x, point_xs[i], y, point_ys[i]) for i in eachindex(points))
    E_tot_mag(x, y, ρs) = hypot(E_tot_x(x, y, ρs), E_tot_y(x, y, ρs))

    fig = Figure()
    ax = Axis(fig[1, 1], aspect=1)
    Φ_tots = [Φ_tot(x, y, ρs) for x in xs, y in ys]
    hm = heatmap!(ax, xs, ys, Φ_tots)
    plt = scatter!(ax, Point2f.(points), marker=[ρ≥0 ? '+' : '-' for ρ in ρs], markersize=40 .* abs.(norm_to_one(ρs)), color=[ρ≥0 ? :red : :blue for ρ in ρs])
    # Colorbar(fig[1, 2], hm, label="Φ_tot")

    E_tot_mags = [E_tot_mag(x, y, ρs) for x in xs, y in ys]
    ax2 = Axis(fig[1, 3], aspect=1)
    hm2 = heatmap!(ax2, xs, ys, E_tot_mags)#, colorrange=(0, 3))
    # Colorbar(fig[1, 4], hm2, label="|E_tot|")

    E_tot_xs = [E_tot_x(x, y, ρs) for x in xs, y in ys]
    ax3 = Axis(fig[2, 1], aspect=1)
    hm3 = heatmap!(ax3, xs, ys, E_tot_xs)#, colorrange=(-3, 3))
    # Colorbar(fig[2, 2], hm3, label="E_tot_x")

    E_tot_ys = [E_tot_y(x, y, ρs) for x in xs, y in ys]
    ax4 = Axis(fig[2, 3], aspect=1)
    hm4 = heatmap!(ax4, xs, ys, E_tot_ys)#, colorrange=(-3, 3))
    # Colorbar(fig[2, 4], hm4, label="E_tot_y")

    # plt = scatter!(Point2f.(points), marker=[ρ≥0 ? '+' : '-' for ρ in ρs], markersize=40 .* abs.(norm_to_one(ρs)), color=[ρ≥0 ? :red : :blue for ρ in ρs])
    
    colsize!(fig.layout, 1, Aspect(1, 1))
    colsize!(fig.layout, 3, Aspect(1, 1))
    resize_to_layout!(fig)
    display(fig)
end
##
let     ##* I have to rotr90 to show a matrix as a matrix
    N = 2
    xs = 1:N
    ys = (N+1):2N
    zs = [1 2; 0 1]|>rotr90
    heatmap(xs, ys, zs)
end

let     ##* With no rot, x is inner index, y is outer
    N = 100
    xs = 1:N
    ys = 1:N
    zs = [x + 10y for x in xs, y in ys]
    # zs = [x^2 * 2y for x in xs, y in ys]
    heatmap(xs, ys, zs)
end