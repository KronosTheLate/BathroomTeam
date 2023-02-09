
#JUST TESTING GIT!!!  -ANDI








using Plots

x=0:0.01:5
y=x.^2
# @show(y)
plot(x,y)

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

plot(t,y[:,1], lw = 5, c=:blue, labels="Mice")
plot!(t,y[:,2], lw = 5, c=:red, labels="Cats", legend=:topleft)
plot!(title = "Lotka Volterra diff.eq.", label = ["Mice" "Cats"])
xlabel!("Time [e.g. years]")
ylabel!("Population []")
plot!(xlims=(0,t_tot),xticks=0:1:t_tot)
plot!(ylims=(0,60))
