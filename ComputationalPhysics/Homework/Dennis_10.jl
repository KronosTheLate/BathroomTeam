if occursin("dennishb", homedir())
    using Pkg
    Pkg.activate("ComputationalPhysics", shared=true)
end
using Trapz  # For numerical integration
using LinearAlgebra
using GLMakie; Makie.inline!(true)
##

function mesh_func(N::Integer)
    xspan = (0, 1)
    mesh = range(xspan[1], xspan[2], N+2)
    return mesh
end
mesh_func(100)

function M_func(mesh::AbstractVector)
    internal_inds = eachindex(mesh)[begin+1:end-1]

    M_diagonal = 1/3 .* [mesh[i+1] - mesh[i-1] for i in internal_inds]
    M_subdiagonal = 1/6 .* [mesh[i+1]-mesh[i] for i in internal_inds[begin+1:end]]
    M_superdiagonal = 1/6 .* [mesh[i+1]-mesh[i] for i in internal_inds[begin:end-1]]
    
    if all(M_subdiagonal .≈ M_superdiagonal)
        return SymTridiagonal(M_diagonal, M_subdiagonal)
    else
        return Tridiagonal(M_subdiagonal, M_diagonal, M_superdiagonal)
    end
end
M_func(mesh_func(100))

function S_func(mesh::AbstractVector, v=x->1)
    internal_inds = eachindex(mesh)[begin+1:end-1]

    #¤ "Assume explicitly that v(x) is constant on each element"
    for i in eachindex(mesh)[begin:end-1]
        if v(mesh[i]+eps()) != v(mesh[i+1]-eps())
            @warn "v(mesh[i]+eps()) != v(mesh[i+1]-eps()) for i=$i. This means that v(x) is not constant on each element."
        end
    end
    S_diagonal = [v(mesh[i]+eps())^2/(mesh[i+1] - mesh[i]) + v(mesh[i]-eps())^2/(mesh[i] - mesh[i-1]) for i in internal_inds]
    S_subdiagonal = [-v(mesh[i]-eps())^2/(mesh[i]-mesh[i-1]) for i in internal_inds[begin+1:end]]
    S_superdiagonal = [-v(mesh[i]+eps())^2/(mesh[i+1]-mesh[i]) for i in internal_inds[begin:end-1]]
    
    # S_diagonal = [v(mesh[i])^2/(mesh[i+1] - mesh[i]) + v(mesh[i])^2/(mesh[i] - mesh[i-1]) for i in internal_inds]
    # S_subdiagonal = [-v(mesh[i])^2/(mesh[i]-mesh[i-1]) for i in internal_inds[begin+1:end]]
    # S_superdiagonal = [-v(mesh[i])^2/(mesh[i+1]-mesh[i]) for i in internal_inds[begin:end-1]]
    
    if all(S_subdiagonal .≈ S_superdiagonal)  # ≈ because in creating the mesh, points can differ slightly due to which numbers are representable by floats
        return SymTridiagonal(S_diagonal, S_subdiagonal)
    else
        return Tridiagonal(S_subdiagonal, S_diagonal, S_superdiagonal)
    end
end
S_func(mesh_func(100))

"""
    Find the field profiles of the eigenmodes of the 
    PDF as posed in task a-c
"""
function main1(N)
    mesh = mesh_func(N)
    M, S = M_func(mesh), S_func(mesh)
    vals, vecs = eigen(S, M)
    vecs = eachcol(vecs)
    return (;M, S, vals, vecs)
end

function norm_to_one(v::AbstractVector)
    most_extreme_val = v[findmax(abs, v)[2]]
    v ./= most_extreme_val
end

let N=10
    scatterlines(main1(N).vecs[1]|>norm_to_one, label="Numerical")
    scatterlines!(sin.(mesh_func(N)[begin+1:end-1]*π), label="Analytical")
    axislegend()
    display(current_figure())
end

if true
    println("N\tResiduals")
    Ns = [10, 20, 50, 100, 200, 500, 1000][begin:end-1]
    # Ns = sort([Ns; Ns.+1])
    # Ns = [round(Int, 10 * 1.2^i) for i in 0:20]
    total_errors_vec = Float64[]
    total_errors_val = Float64[]
    for N in Ns
        fig = Figure()
        ax = Axis(fig[1, 1], title="N = $N")
        mesh_interior = mesh_func(N)[begin+1:end-1]
        sol = main1(N)
        numerical = sol.vecs[1]
        numerical .= norm_to_one(numerical)
        # most_extreme_val = numerical[findmax(abs, numerical)[2]]
        # numerical ./= most_extreme_val
        reference = sin.(mesh_interior .* π)
        residuals = numerical - reference
        # total_error_vec = norm(residuals, Inf)
        total_error_vec = √trapz(mesh_interior, residuals.^2)

        # log10.(numerical)
        # log10.(reference)
        # log10.(residuals)  #¤ This errors

        scatterlines!(ax, mesh_interior, numerical, label="Numerical")
        scatterlines!(ax, mesh_interior, reference, label="Reference")
        scatterlines!(ax, mesh_interior, residuals .+ 1e-10, label="Residuals + 10⁻¹⁰", color=:red)
        axislegend(position=:rb)
        
        push!(total_errors_vec, total_error_vec)
        push!(total_errors_val, sol.vals[1]-π^2)
        println(N, "\t", norm(residuals))

        try 
            # ylims!(ax, 1e-8, 1)
            ax.yscale = log10
            display(fig)
        catch e
            @show e
            @show extrema(residuals)
        end
        # ax.xscale = log10
    end
    scatterlines(Ns, total_errors_vec, label="Eigenvector error",
        axis=(xlabel="N", xscale=log10, yscale=log10)
    )
    scatterlines!(Ns, total_errors_val, label="Eigenvalue error")
    axislegend(position=:lb)
    display(current_figure())
end

##¤ d)

"""
    Find the field profiles of the eigenmodes of the 
    PDF as posed in task d).
    N is the number of elements at each side of the `boundry_node`
"""
function main2(N)
    v²(x) = x < √0.5 ? 1 : 0.1
    v(x) = √v²(x)
    bndry_node = √0.5
    mesh_part1 = range(0, bndry_node, N+1)
    mesh_part2 = range(bndry_node, 1, N+1)
    mesh = [mesh_part1 ; mesh_part2[begin+1:end]]
    M, S = M_func(mesh), S_func(mesh, v)
    vals, vecs = eigen(S, M)
    vecs = eachcol(vecs)
    return (;bndry_node, mesh, M, S, vals, vecs)
end

let 
    N = 21
    v²(x) = x < √0.5 ? 1 : 0.1
    v(x) = √v²(x)
    bndry_node = √0.5
    mesh_part1 = range(0, bndry_node, N+1)
    mesh_part2 = range(bndry_node, 1, N+1)
    mesh = [mesh_part1 ; mesh_part2[begin+1:end]]
    vscodedisplay(S_func(mesh, v))
end

let N = 201
    sol = main2(N)
    @show count(<(√0.5), sol.mesh)
    @show count(>(√0.5), sol.mesh)
    mesh_interior = sol.mesh[begin+1:end-1]
    sol.vecs[1]
    print("ω = ")
    @show √sol.vals[1]
    fig = Figure()
    ax = Axis(fig[1, 1], title="$N elements on each side, $(length(mesh_interior)) interior points total")
    numerical = sol.vecs[1]
    most_extreme_val = numerical[findmax(abs, numerical)[2]]
    numerical ./= most_extreme_val
    scatterlines!(ax, sol.mesh, [0; numerical; 0], label="Numerical")
    vlines!([√0.5])
    display(fig)
    # most_extreme_val = numerical[findmax(abs, numerical)[2]]
    # numerical ./= most_extreme_val
    # reference = sin.(mesh_interior .* π)
    # residuals = numerical - reference
    # total_error_vec = √trapz(mesh_interior, residuals.^2)
end