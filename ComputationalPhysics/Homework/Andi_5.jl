using LinearAlgebra

# Only fitted to Keplar problem, since knowledge of underlying variables seems needed.
function leap_frog(f,x₀,t)
    x=copy(x₀)
    h=step(t)
    xtot=zeros(size(x,1),size(t,1))
    xtot[:,1]=x
    for i in 2:size(t,1)
        dx=f(x)
        x[1:2]+=dx[1:2]*h # Could interchange position and momentum update. 1:2->3:4
        dx=f(x)
        x[3:4]+=dx[3:4]*h # Could interchange position and momentum update. 3:4->1:2
        xtot[:,i]=x
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

    return hcat(drx, dry, drz, dpx, dpy, dpz)
end


##


t₀=0
t₁=50
Δt=0.01
Nt=Int(round((t₁-t₀)/Δt))
t=range(t₀,t₁,Nt)


N=300
# rx = Vector{N}
rx = Vector{Float64}(undef,N)
# rx[3] = 1
# display(rx)
ry = Vector{Float64}(undef,N)
rz = Vector{Float64}(undef,N)
px = Vector{Float64}(undef,N)
py = Vector{Float64}(undef,N)
pz = Vector{Float64}(undef,N)

states=[1 1 1 1 1 1 ; 2 2 2 2 2 2]
change=d_dt_kepler(states,[0;0;0])
display(change)


cat(reshape(states,2,6,1),reshape(change,2,6,1),dims=3)



