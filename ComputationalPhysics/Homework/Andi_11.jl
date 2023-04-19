
using LinearAlgebra 
using Plots
# plotly() # Interactive plots
gr()


function 𝒢_Φ(x,x′,y,y′)
    return @. -1/(4π)*log((x-x′)^2+(y-y′)^2)
end

function 𝒢_Ex(x,x′,y,y′)
    return -1/(4π)*1/((x-x′)^2+(y-y′)^2) * 2*(x-x′)
end

function 𝒢_Ey(x,x′,y,y′)
    return -1/(4π)*1/((x-x′)^2+(y-y′)^2) * 2*(y-y′)
end

function Eabsᵢₙₜ()
    return error("Not implemented. Maybe not necessary for function")
end

# Note the usage of a and b in several of the below functions.
a=1
b=0
function Φₑₓₜ(x,y)
    return a*x+b*y
end

function Exₑₓₜ(x,y)
    return a
end

function Eyₑₓₜ(x,y)
    return b
end


# NN=101

xs = range(-2,2,NN)

# xsm = vcat([repeat(xs[i,:]', inner=(r[i],1)) for i in indices(xs,1)]...)


ys = range(-2,2,NN)

x′=1
y′=0

bla = heatmap(xs, ys, 𝒢_Φ(xs,x′,ys,y′),xlims=(-2,2), ylims=(-2,2))



# heatmap(x,y,𝒢_Φ(x,x′,y,y′),xlims=(-2,2), ylims=(-2,2))



