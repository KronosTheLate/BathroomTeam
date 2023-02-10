
#JUST TESTING GIT!!!  -ANDI






function propagationM(n,D,k)
    P=[exp(im*n*D*k) 0;
       0 exp(im*n*D*k)]
    return P
end

function interfaceM(nᵢ₋₁,nᵢ)
    nr=nᵢ/nᵢ₋₁
    I=1/2*[1+nr   1-nr;
           1-nr   1+nr]
    return I
end


#@show(propagationM(1.5,1,1))
#@show(interfaceM(1,2))



##

using Plots
plotly()

x=0:0.01:5
y=x.^2
# @show(y)
Plots.plot(x,y)

## Lotka-Volterra coupled differential equation

function LotkaVolterra(y::Vector,dt;α=1,β=0.1,γ=2,δ=0.05)
    🐁=y[1]
    🐈=y[2]
    d🐁  = (α*🐁 - β*🐁*🐈)*dt
    d🐈 = (-γ*🐈 + δ*🐁*🐈)*dt
    return dy=[d🐁, d🐈]
end



dt=0.001
N=Int.(10/dt) 
t_tot=dt*N
t=0:dt:t_tot-dt
y=zeros(N,2)

y_init=[50,15]

y[1,:]=y_init
for i in 1:(N-1)
    y[i+1,:]=y[i,:]+LotkaVolterra(y[i,:],dt)
end
# @show(y)

Plots.plot(t,y[:,1], lw = 5, c=:blue, labels="Mice")
Plots.plot!(t,y[:,2], lw = 5, c=:red, labels="Cats", legend=:topleft)
Plots.plot!(title = "Lotka Volterra diff.eq.", label = ["Mice" "Cats"])
Plots.xlabel!("Time [e.g. years]")
Plots.ylabel!("Population []")
Plots.plot!(xlims=(0,t_tot),xticks=0:1:t_tot)
Plots.plot!(ylims=(0,60))
