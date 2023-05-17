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
    ax3, hm3 = heatmap(fig[2, 1], xs, ys, 𝒢_Exs, axis=(title="E_x", aspect=1, xlabel="x", ylabel="y"))
    Colorbar(fig[2, 2], hm3)
    ax4, hm4 = heatmap(fig[2, 3], xs, ys, 𝒢_Eys, axis=(title="E_y", aspect=1, xlabel="x", ylabel="y"))
    Colorbar(fig[2, 4], hm4)
    colsize!(fig.layout, 1, Aspect(1, 1))
    colsize!(fig.layout, 3, Aspect(1, 1))
    resize_to_layout!(fig)
    display(fig)
end

##¤ a2)
# Φ_ext_func(p) = [p.a*x + p.b*y for x in getindex.(p.points, 1), y in getindex.(p.points, 2)]
Φ_ext_func(p) = [p.a*x + p.b*y for x in p.xs, y in p.ys]
E_ext_x_func(p) = [-p.a for x in p.xs, y in p.ys]
E_ext_y_func(p) = [-p.b for x in p.xs, y in p.ys]
Φ_ext_atpoints_func(p) = [p.a*x + p.b*y for (x, y) in p.points]
E_ext_x_atpoints_func(p) = [-p.a for (x, y) in p.points]
E_ext_y_atpoints_func(p) = [-p.b for (x, y) in p.points]


#¤ Other funcs
Φ_int_func(ρs, p) = [sum(ρs[i]*𝒢_Φ(x, p.points[i][1], y, p.points[i][2]) for i in eachindex(p.points)) for x in p.xs, y in p.ys]
E_int_x_func(ρs, p) = [sum(ρs[i]*𝒢_Ex(x, p.points[i][1], y, p.points[i][2]) for i in eachindex(p.points)) for x in p.xs, y in p.ys]
E_int_y_func(ρs, p) = [sum(ρs[i]*𝒢_Ey(x, p.points[i][1], y, p.points[i][2]) for i in eachindex(p.points)) for x in p.xs, y in p.ys]

Φ_tot_func(ρs, p) = Φ_ext_func(p) .+ Φ_int_func(ρs, p)
E_tot_x_func(ρs, p) = E_ext_x_func(p) + E_int_x_func(ρs, p)
E_tot_y_func(ρs, p) = E_ext_y_func(p) + E_int_y_func(ρs, p)
E_tot_mag_func(ρs, p) = hypot.(E_tot_x_func(ρs, p), E_tot_y_func(ρs, p))

let     
    #¤ b)
    p = parameters(N_mesh=100, N_points=100, b=0)
    (;N_mesh, N_points, xs, ys, points) = p
    point_xs = getindex.(p.points, 1)
    point_ys = getindex.(p.points, 2)
    b⃗ = -Φ_ext_atpoints_func(p)
    M = [i==j ? 0.0 : 𝒢_Φ(point_xs[i], point_xs[j], point_ys[i], point_ys[j]) for i in eachindex(points), j in eachindex(points)]# .|> ComplexF64
    M = issymmetric(M) ? Symmetric(M) : M

    ρs = M \ b⃗

    Φ_int = Φ_int_func(ρs, p)
    Φ_ext = Φ_ext_func(p)
    Φ_tot = Φ_tot_func(ρs, p)

    E_tot_mag = E_tot_mag_func(ρs, p)
    E_tot_x = E_tot_x_func(ρs, p)
    E_tot_y = E_tot_y_func(ρs, p)

    let
        fig = Figure()
        Label(fig[0, 1:4], rich("Potential plots", font=:bold, fontsize=30))
        
        ax = Axis(fig[1, 1], aspect=1)
        hm = heatmap!(ax, xs, ys, Φ_tot)
        # return Point2f.(p.points)
        plt = scatter!(ax, Point2f.(p.points), marker=[ρ≥0 ? '+' : '-' for ρ in ρs], markersize=40 .* abs.(normalize(ρs, Inf)), color=[ρ≥0 ? :red : :blue for ρ in ρs])
        Colorbar(fig[1, 2], hm, label="Φ_tot")

        ax2 = Axis(fig[1, 3], aspect=1)
        # hm2 = heatmap!(ax2, xs, ys, E_tot_mag, colorrange=(0, 3))
        # Colorbar(fig[1, 4], hm2, label="|E_tot|")
        hm2 = heatmap!(ax2, xs, ys, Φ_tot.|>abs.|>log10)#, colorrange=(0, 3))
        Colorbar(fig[1, 4], hm2, label="log(|Φ_tot|)")

        ax3 = Axis(fig[2, 1], aspect=1)
        # hm3 = heatmap!(ax3, xs, ys, E_tot_x, colorrange=(-3, 3))
        # Colorbar(fig[2, 2], hm3, label="E_tot_x")
        hm3 = heatmap!(ax3, xs, ys, Φ_ext, colorrange=(-3, 3))
        Colorbar(fig[2, 2], hm3, label="Φ_ext")

        ax4 = Axis(fig[2, 3], aspect=1)
        # hm4 = heatmap!(ax4, xs, ys, E_tot_y, colorrange=(-3, 3))
        # Colorbar(fig[2, 4], hm4, label="E_tot_y")
        hm4 = heatmap!(ax4, xs, ys, Φ_int, colorrange=(-3, 3))
        Colorbar(fig[2, 4], hm4, label="Φ_int")
        colsize!(fig.layout, 1, Aspect(1, 1))
        colsize!(fig.layout, 3, Aspect(1, 1))
        resize_to_layout!(fig)
        display(fig)
    end

    let
        fig = Figure()
        Label(fig[0, 1:4], rich("E-field plots", font=:bold, fontsize=30))
        
        ax = Axis(fig[1, 1], aspect=1)
        hm = heatmap!(ax, xs, ys, E_tot_mag, colorrange=(0, 3))
        plt = scatter!(ax, Point2f.(points), marker=[ρ≥0 ? '+' : '-' for ρ in ρs], markersize=40 .* abs.(normalize(ρs, Inf)), color=[ρ≥0 ? :red : :blue for ρ in ρs])
        Colorbar(fig[1, 2], hm, label="|E_tot|")

        ax2 = Axis(fig[1, 3], aspect=1)
        hm2 = heatmap!(ax2, xs, ys, E_tot_mag.|>log10)
        Colorbar(fig[1, 4], hm2, label="log(|E_tot|)")

        ax3 = Axis(fig[2, 1], aspect=1)
        hm3 = heatmap!(ax3, xs, ys, E_tot_x, colorrange=(-3, 3))
        Colorbar(fig[2, 2], hm3, label="E_tot_x")

        ax4 = Axis(fig[2, 3], aspect=1)
        hm4 = heatmap!(ax4, xs, ys, E_tot_y, colorrange=(-3, 3))
        Colorbar(fig[2, 4], hm4, label="E_tot_y")

        colsize!(fig.layout, 1, Aspect(1, 1))
        colsize!(fig.layout, 3, Aspect(1, 1))
        resize_to_layout!(fig)
        display(fig)
    end
end

##¤ Convergence