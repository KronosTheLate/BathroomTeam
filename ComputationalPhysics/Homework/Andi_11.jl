
using LinearAlgebra 
using Plots
# plotly() # Interactive plots
gr()


function 𝒢_Φ(x,x′,y,y′)
    return @. -1/(4π)*log((x-x′)^2+(y-y′)^2)
end

function 𝒢_Ex(x,x′,y,y′)
    return @. 1/(4π)*1/((x-x′)^2+(y-y′)^2) * 2*(x-x′)
end

function 𝒢_Ey(x,x′,y,y′)
    return @. 1/(4π)*1/((x-x′)^2+(y-y′)^2) * 2*(y-y′)
end

function Eabsᵢₙₜ(x,x′,y,y′)
    return @. √(𝒢_Ex(x,x′,y,y′)^2+𝒢_Ey(x,x′,y,y′)^2)
    # throw(error())
    # error("Not implemented. Maybe not necessary for function")
end

# Note the usage of a and b in several of the below functions.
a1=1
b1=0
function Φₑₓₜ(x,y)
    return @. a1*x+b1*y
end

function Exₑₓₜ(x,y)
    return @. -a1
end

function Eyₑₓₜ(x,y)
    return @. -b1
end


NN=200

xs = range(-2,2,NN)
ys = range(-2,2,NN)

x′=1
y′=0
𝒢s = [𝒢_Φ(x, x′, y, y′) for x in xs, y in ys]
𝒢s_plot = heatmap(xs, ys, reverse(rotl90(𝒢s), dims = 1))
display(𝒢s_plot)

Exs = [𝒢_Ex(x,x′,y,y′) for x in xs, y in ys]
Exs_plot = heatmap(xs, ys, reverse(rotl90(Exs), dims = 1), clim=(-5,5)) #clim=(-5,5)
display(Exs_plot)

Eys = [𝒢_Ey(x,x′,y,y′) for x in xs, y in ys]
Eys_plot = heatmap(xs, ys, reverse(rotl90(Eys), dims = 1), clim=(-5,5)) #clim=(-5,5)
# scatter!(Eys_plot, [1], [0.5], markersize=2, markercolor=:white)
        # reverse(rotl90(Eys), dims = 1)
display(Eys_plot)

Eabss = [Eabsᵢₙₜ(x,x′,y,y′) for x in xs, y in ys]
Eabss_plot = heatmap(xs, ys, reverse(rotl90(Eabss), dims = 1), clim=(0,5)) #, clim=(0,5)
display(Eabss_plot)

Φs = [Φₑₓₜ(x,y) for x in xs, y in ys]
Φs_plot = heatmap(xs, ys, reverse(rotl90(Φs), dims = 1))
display(Φs_plot)


# heatmap([1 2 3;4 5 6;7 8 9]) # To see how to rotate/orient matrix

##

N=100
θ=range(0,2π,N+1)
θ=θ[1:end-1]

rxs = cos.(θ)
rys = sin.(θ)
rs = hcat(rxs,rys)

charge_pos_plot=scatter(rs[:,1],rs[:,2],aspect_ratio = :equal, markershape=:circle, label="Charge positions")
display(charge_pos_plot)


b = -Φₑₓₜ(rxs,rys)

M =  [i==j ? 0.0 : 𝒢_Φ(rxs[i], rxs[j], rys[i], rys[j]) for i in eachindex(rxs), j in eachindex(rys)] 
# M =  [i==j ? 0.0 : 1*i+0*j for i in eachindex(rxs), j in eachindex(rys)] 
    # Note the indexes. Just have to loop over 
# display(M)

M = issymmetric(M) ? Symmetric(M) : M

ρ = M\b
charge_pos_plot=scatter(rs[:,1],rs[:,2],aspect_ratio = :equal, markershape=:circle, label="Charge positions", zcolor=ρ)
display(charge_pos_plot)



Φₜₒₜ = Φs
for i in 1:size(rxs,1)
    # print(i)
    x′ = rxs[i]
    y′ = rys[i]
    Φₜₒₜ += ρ[i].*[𝒢_Φ(x, x′, y, y′) for x in xs, y in ys]
end

# NOTE the circle marker color.
Φₜₒₜ_plot = heatmap(xs, ys, reverse(rotl90(Φₜₒₜ), dims = 1), aspect_ratio=:equal)
scatter!(Φₜₒₜ_plot, rs[:,1],rs[:,2], markershape=:circle, markersize=1.5, label="", zcolor=ρ.*(N/15)) #, markercolor = :white


##

Exs_tot = [Exₑₓₜ(x, y) for x in xs, y in ys]
for i in 1:size(rxs,1)
    # print(i)
    x′ = rxs[i]
    y′ = rys[i]
    Exs_tot += ρ[i].*[𝒢_Ex(x, x′, y, y′) for x in xs, y in ys]
end

Exs_tot_plot = heatmap(xs, ys, reverse(rotl90(Exs_tot), dims = 1), aspect_ratio=:equal)
scatter!(Exs_tot_plot, rs[:,1],rs[:,2], markershape=:circle, markersize=1.5, markercolor = :white , label="")#, zcolor=ρ.*(N/15))
display(Exs_tot_plot)


Eys_tot = [Eyₑₓₜ(x, y) for x in xs, y in ys]
for i in 1:size(rxs,1)
    # print(i)
    x′ = rxs[i]
    y′ = rys[i]
    Eys_tot += ρ[i].*[𝒢_Ey(x, x′, y, y′) for x in xs, y in ys]
end

# Note the clim, since one point explode in value
Eys_tot_plot = heatmap(xs, ys, reverse(rotl90(Eys_tot), dims = 1), aspect_ratio=:equal) #clim=(-1.5,1.5)
scatter!(Eys_tot_plot, rs[:,1],rs[:,2], markershape=:circle, markersize=1.5, markercolor = :white , label="")#, zcolor=ρ.*(N/15))
display(Eys_tot_plot)


E_tot = @. sqrt(Exs_tot^2+Eys_tot^2)

E_tot_plot = heatmap(xs, ys, reverse(rotl90(E_tot), dims = 1), aspect_ratio=:equal,clim=(0,3)) #clim=(-1.5,1.5)
scatter!(E_tot_plot, rs[:,1],rs[:,2], markershape=:circle, markersize=1.5, markercolor = :white , label="")#, zcolor=ρ.*(N/15))
display(E_tot_plot)


