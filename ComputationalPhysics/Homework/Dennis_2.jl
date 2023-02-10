using LinearAlgebra  # To get identity matrix of arbitrary size
using GLMakie; Makie.inline!(true)  # Plotting package. inline to use plot pane. Say false to get pop-outs

##

"""
Return propagation transformation matrix
"""
P(n, D, k₀) = [cis(n*D*k₀) 0; 0 cis(-n*D*k₀)]

"""
Return interface transformation matrix
"""
I(n₁, n₂) = 1/2 * [
    1+n₂/n₁     1-n₂/n₁
    1-n₂/n₁     1+n₂/n₁
]

function transfer_matrix(n⃗, D⃗, k₀)
    # Assume first and last transformations are transmissions, with propagation between each.
    # This means that you propegate D⃗[1] in medium of n = n⃗[2], as no propagation occurs in n⃗[1]

    # Need 2 values of n for initial and final medium, and both n & D for mediums in between.
    (length(n⃗) == length(D⃗)+2)  ||  throw(DimensionMismatch("Expected `x+2` refractive indices for `x` distances.\nGot $(length(n⃗)) refractive indices and $(length(D⃗)) distances."))
    tm = I(n⃗[1], n⃗[2])
    for i in eachindex(D⃗)
        tm = P(n⃗[i+1], D⃗[i], k₀) * tm  # propagation
        tm = I(n⃗[i+1], n⃗[i+2])   * tm    # interface
    end
    tm
end

function transmittance(args...)
    tm_ = transfer_matrix(args...)
    # println("tm_ = ")  # printing for troubleshooting...
    # display(tm_)       # values are complex and changing
    return abs2(tm_[1, 1] - tm_[1, 2]*tm_[2, 1]/tm_[2, 2])
end

function reflectance(args...)
    tm_ = transfer_matrix(args...)
    return abs2(tm_[2, 1]/tm_[2, 2])
end

##¤ d) part 1
with_theme(resolution=(1920÷2, 1080÷2.2)) do
    for n_interfaces in (11, 21, 41)
        n_propagations = n_interfaces-2               # Just being very explicit for myself
        n⃗s = [iseven(i) ? 1 : 2 for i in 1:n_interfaces]
        D⃗s = ones(n_propagations)
        k₀s = range(0, 3, 1000)

        transmittances = [transmittance(n⃗s, D⃗s, k₀) for k₀ in k₀s]
        reflectances = [reflectance(n⃗s, D⃗s, k₀) for k₀ in k₀s]
        E_lost = 1 .- transmittances .- reflectances

        fig, ax, plt = lines(k₀s, transmittances, label="Transmittance")
        lines!(k₀s, reflectances, label="Reflectance")
        axislegend(position=(1, 0.5))
        lines(fig[2, 1], k₀s, E_lost, label=L"E_{lost}")
        axislegend(position=(1, 0))
        ax.xlabel = "k₀"
        ax.ylabel = "Value"
        ax.title = "$n_interfaces interfaces and $n_propagations propagations"
        fig |> display
        # transmittances|>display
    end
end

##¤ d) part 2
with_theme(resolution=(1920÷2, 1080÷2.2)) do
    let n_interfaces = 11
        n_propagations = n_interfaces-2               # Just being very explicit for myself
        n⃗s = [iseven(i) ? 1 : 2 for i in 1:n_interfaces]
        n⃗s[end÷2] += 3
        n⃗s[end÷2+3] += 3
        D⃗s = ones(Int64, n_propagations)
        D⃗s[2] += 17
        k₀s = range(0, 3, 1000)

        transmittances = [transmittance(n⃗s, D⃗s, k₀) for k₀ in k₀s]
        reflectances = [reflectance(n⃗s, D⃗s, k₀) for k₀ in k₀s]
        E_lost = 1 .- transmittances .- reflectances

        fig, ax, plt = lines(k₀s, transmittances, label="Transmittance")
        lines!(k₀s, reflectances, label="Reflectance")
        axislegend(position=(1, 0.5))
        lines(fig[2, 1], k₀s, E_lost, label=L"E_{lost}")
        axislegend(position=(1, 0))
        Label(fig[end+1, :], "Indices of refraction:   "* join(n⃗s, "   "), tellwidth=false, textsize=24)
        Label(fig[end+1, :], rpad("Thicknesses:", 27)* join(D⃗s, "   "), tellwidth=false, textsize=24)
        ax.xlabel = "k₀"
        ax.ylabel = "Value"
        ax.title = "$n_interfaces interfaces and $n_propagations propagations"
        fig |> display
        # transmittances|>display
    end
end

##¤ e)

with_theme(resolution=(1920÷2, 1080÷2.2)) do
    for n_central in (2+5im, 2+10im, 2+15im)
        n_interfaces = 11
        n_propagations = n_interfaces-2               # Just being very explicit for myself
        n⃗s = Any[iseven(i) ? 1 : 2 for i in 1:n_interfaces]
        n⃗s[end÷2] = n_central
        D⃗s = ones(Int64, n_propagations)
        k₀s = range(0, 3, 1000)

        transmittances = [transmittance(n⃗s, D⃗s, k₀) for k₀ in k₀s]
        reflectances = [reflectance(n⃗s, D⃗s, k₀) for k₀ in k₀s]
        E_lost = 1 .- transmittances .- reflectances

        fig, ax, plt = lines(k₀s, transmittances, label="Transmittance")
        lines!(k₀s, reflectances, label="Reflectance")
        axislegend(position=(1, 0.5))
        lines(fig[2, 1], k₀s, E_lost, label=L"E_{lost}")
        axislegend(position=(1, 0))
        Label(fig[end+1, :], "Indices of refraction:   "* join(n⃗s, "   "), tellwidth=false, textsize=24)
        Label(fig[end+1, :], rpad("Thicknesses:", 27)* join(D⃗s, "   "), tellwidth=false, textsize=24)
        ax.xlabel = "k₀"
        ax.ylabel = "Value"
        ax.title = "$n_interfaces interfaces and $n_propagations propagations"
        fig |> display
        # transmittances|>display
    end
end