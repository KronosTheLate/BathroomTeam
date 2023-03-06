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
function d_dt_kepler(rp)
    #rp of format rp=[x,y,p_x,p_y], and similar for drp
    drp=[0.0; 0.0; 0.0; 0.0] #initialize to index
    # See diff.eq. in problem sheet
    drp[1]=rp[3]
    drp[2]=rp[4]
    drp[3]=-1/norm(rp[1:2])^3*rp[1]
    drp[4]=-1/norm(rp[1:2])^3*rp[2]
    return drp
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
rx[3] = 1
# display(rx)
ry = Vector{Float64}(undef,N)
rz = Vector{Float64}(undef,N)
px = Vector{Float64}(undef,N)
py = Vector{Float64}(undef,N)
pz = Vector{Float64}(undef,N)


