using LinearAlgebra
using Random # random number in "rand"
using Distributions # uniform distribution with "Uniform()"
using Plots
# plotly() # Interactive plots
gr()
# using PyPlot


# Only fitted to Keplar problem, since knowledge of underlying variables seems needed.
function leap_frog(f,x₀,t,R)
    x=copy(x₀)
    h=step(t)
    xtot=zeros(size(x,1),size(x,2),size(t,1))
    xtot[:,:,1]=x
    for i in 2:size(t,1)
        if t[i]>25 R=[1;0;0] end #Quick and dirty task 1b for particle box moving
        dx=f(x,R)
        # display(x)
        # display(dx)
        x[:,1:3]+=dx[:,1:3]*h # Could interchange position and momentum update. 1:3->4:6
        # display(x)
        dx=f(x,R)
        x[:,4:6]+=dx[:,4:6]*h # Could interchange position and momentum update. 4:6->1:3
        xtot[:,:,i]=x
    end
    return xtot
end
leap_frog(args)=leap_frog(args...)

# Differential equations to describe Kepler problem of one object orbiting 0,0 , like heavy sun.
function d_dt_kepler(states,R)
    N=size(states,1)

    rx=states[:,1]
    ry=states[:,2]
    rz=states[:,3]
    px=states[:,4]
    py=states[:,5]
    pz=states[:,6]

    drx = Vector{Float64}(undef,N)
    dry = Vector{Float64}(undef,N)
    drz = Vector{Float64}(undef,N)
    dpx = Vector{Float64}(undef,N)
    dpy = Vector{Float64}(undef,N)
    dpz = Vector{Float64}(undef,N)
    # dp  = Vector{Vector{Float64}}(undef,N)

    m=1
    D=1
    # See diff.eq. in problem sheet

    for i in 1:N
        drx[i]=1/m*px[i]
        dry[i]=1/m*py[i]
        drz[i]=1/m*pz[i]

        dp=-D*norm([rx[i]; ry[i]; rz[i]]-R)^2 .*([rx[i]; ry[i]; rz[i]]-R)
        
        # dp_temp=dp[i]
        dpx[i]=dp[1]
        dpy[i]=dp[2]
        dpz[i]=dp[3]
    end

    if true # Repulsive forces
        #Implementation of repulsive force between particles. Can be toggled on and off
        F=zeros(N,3)
        α=1/N
        for j in 1:N #Note could cut computation time down, since every computation is calculated twice. 
                     #E.g. j=2 and jj=3 is similar but with *-1 as j=3 and jj=2
            rj=[rx[j];ry[j];rz[j]]
            for jj in 1:N
                if jj != j        
                    rjj=[rx[jj];ry[jj];rz[jj]]
                    F[j,:]+=α*(rj-rjj) / norm(rj-rjj)^3
                end
            end
            dpx[j]+=F[j,1]
            dpy[j]+=F[j,2]
            dpz[j]+=F[j,3]
        end
    end

    return hcat(drx, dry, drz, dpx, dpy, dpz)
end


##


t₀=0
t₁=50
Δt=0.01
Nt=Int(round((t₁-t₀)/Δt))
t=range(t₀,t₁,Nt)


N=30



rng = MersenneTwister(1234)
# rx = Vector{N}
rx = rand(rng, Uniform(-0.5,0.5),N) #Vector{Float64}(undef,N)
# rand!(rng, Uniform(-0.5,0.5), rx)
# rx[3] = 1
# display(rx)
ry = rand(rng, Uniform(-0.5,0.5),N) #Vector{Float64}(undef,N)
rz = rand(rng, Uniform(-0.5,0.5),N) #Vector{Float64}(undef,N)
px = rand(rng, Uniform(-0.5,0.5),N) #Vector{Float64}(undef,N)
py = rand(rng, Uniform(-0.5,0.5),N) #Vector{Float64}(undef,N)
pz = rand(rng, Uniform(-0.5,0.5),N) #Vector{Float64}(undef,N)

states₀=hcat(rx,ry,rz,px,py,pz)

R=[0; 0; 0]

states_total=leap_frog(d_dt_kepler,states₀,t,R)


rx=states_total[:,1,:]
ry=states_total[:,2,:]
rz=states_total[:,3,:]

##

t_range=range(t₀,t₁,Nt)
#https://discourse.julialang.org/t/3d-animation-plot-of-2-ode-solutions/83408
nframes = 200 # Change later for large Nt
# anim = @animate for t in LinRange(first(t), last(t), nframes)
anim = @animate for j in LinRange(1, size(rx,2), nframes)
    i=Int(round(j))
    scatter(rx[:,i],ry[:,i],rz[:,i], legend=false)#, color=3, lab = "time=$(round(t, digits = 2))")
    # plot!(sol2, vars = (1,2,3), tspan = (0.0, t), lab = "Solution 2")
    
    #current time marker
    # scatter!(rx[:,i],ry[:,i],rz[:,i], vars = (1,2,3), tspan = (t, t), lab = nothing, color = 1)
    # scatter!(sol2, vars = (1,2,3), tspan = (t, t), lab = nothing, color = 2)
    
    # axis setting... need to be last b/c DiffEq recipes will overwrite.
    # plot!(xaxis = ("x" ,(-30, 30)), yaxis = ("y", (-30,30)), zaxis=("z", (0, 60)),
        # title = "t = $(round(t, digits = 2))")
    view=3
    plot!(xaxis = ("x" ,(-view, view)), yaxis = ("y", (-view,view)), zaxis=("z", (-view, view)),title = "N=$N, t = $(round(t[i], digits = 2))")
end

gif(anim, "C:\\Users\\andir\\Downloads\\Andi_5_subtask_a.gif"; fps = 20)



