using LinearAlgebra  # To get identity matrix of arbitrary size
using GLMakie; Makie.inline!(true)  # Plotting package. inline to use plot pane. Say false to get pop-outs

##

"""
Return propagation transformation matrix
"""
P(n, D, k₀) = [cis(-n*D*k₀) 0; 0 cis(-n*D*k₀)]

#Andi's def:
P(n, D, k₀) = [cis(n*D*k₀) 0; 0 -cis(n*D*k₀)]

"""
Return interface transformation matrix
"""
I(n₁, n₂) = 1/2 * [
    1+n₂/n₁     1-n₂/n₁
    1-n₂/n₁     1+n₂/n₁
]

"""
Given a vector of refractive indexes n⃗ and 
distances D⃗, return transfer matrix.
"""
function blafunc(n⃗, D⃗, k₀)
    # Assume first and last transformations are transmissions, with propagation between each.
    # This means that you propegate D⃗[1] in medium of n = n⃗[2], as no propagation occurs in n⃗[1]

    # Need 2 values of n for initial and final medium, and both n & D for mediums in between.
    (length(n⃗) == length(D⃗)+2)  ||  throw(ArgumentError("Expected `x+2` refractive indices for `x` distances.\nGot $(length(n⃗)) refractive indices and $(length(D⃗)) distances."))
    transfer_matrix = I(n⃗[1], n⃗[2])
    for i in eachindex(D⃗)
        transfer_matrix = P(n⃗[i+1], D⃗[i], k₀) * transfer_matrix  # propagation
        transfer_matrix = I(n⃗[i+1], n⃗[i+2])   * transfer_matrix    # interface
    end
    println("Displaying transfer_matrix")
    display(transfer_matrix)
    @show typeof(transfer_matrix)
    transfer_matrix
end
blafunc([1, 2, 3], [4], 5)


function T()
    my_matrix = rand(2, 2)
    println("Displaying my_matrix")
    display(my_matrix)
    @show typeof(my_matrix)
    return my_matrix
end
𝚃()
𝚃([1, 2, 3], [4], 5)
P(1, 1, 1)
"""
Transmittance. Typed by \\ttT<tab>
"""
function 𝚃(args...)
    T_ = T(args...)
    # println("T_ = ")  # printing for troubleshooting...
    # display(T_)       # values are complex and changing
    return abs2(T_[1, 1] - T_[1, 2]*T_[2, 1]/T_[2, 2])
end

"""
Reflectance. Typed by \\ttR<tab>
"""
function 𝚁(args...)
    T_ = T(args...)
    return abs2(T_[2, 1]/T_[2, 2])
end

##

let # using let block for namespace hygiene
    n_interfaces = 10                   # or 20, or 40      
    n_propagations = n_interfaces-2               # Just being very explicit for myself
    n⃗s = [iseven(i) ? 1 : 2 for i in 1:n_interfaces]
    D⃗s = ones(n_propagations)
    k₀s = range(0, 3, 1000)

    𝚃s = [𝚃(n⃗s, D⃗s, k₀) for k₀ in k₀s]
    𝚁s = [𝚁(n⃗s, D⃗s, k₀) for k₀ in k₀s]
    E_lost = 1 .- 𝚃s .- 𝚁s

    fig, ax, plt = scatterlines(k₀s, 𝚃s, label="𝚃")
    # scatterlines!(k₀s, 𝚁s, label="𝚁")
    # scatterlines!(k₀s, E_lost, label=L"E_{lost}")
    Legend(fig[1, 2], ax)
    ax.xlabel = "k₀"
    ax.ylabel = "Value"
    fig #|> display
    𝚃s|>display
end


let # Verifying that abs2 of a complex number gives same as norm squared
    x = 17*cis(deg2rad(77))
    norm(x)^2, abs2(x)
end