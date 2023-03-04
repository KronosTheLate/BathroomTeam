using LinearAlgebra 
using LaTeXStrings
using GLMakie; Makie.inline!(true)

# More general Euler method. Used for task 2
function euler(f,x₀,t)
    x=x₀
    h=step(t)
    xtot=zeros(size(x,1),size(t,1))
    xtot[:,1]=x
    for i in 2:size(t,1)
        dx=f(x)
        x+=dx*h
        xtot[:,i]=x
    end
    
    return xtot
end
euler(args)=euler(args...)

function plot_orbit2(xtable,t) # Credit Dennis. Changed to fit my implementation -Andi
    # begin  # Plotting components
        r₀=xtable[1,1]
        p₀=xtable[4,1]
        N=size(t,1)
        r_xs=xtable[1,:]
        r_ys=xtable[2,:]
        p_xs=xtable[3,:]
        p_ys=xtable[4,:]
        fig = Figure()
        ax1 = Axis(fig[1, 1], title="r₀ = $r₀\np⃗₀ = $p₀\nN = $N", xlabel="Timepoint", ylabel="X component",)
        scatterlines!(t, r_xs, label=L"\vec{r}_x", color = Cycled(1)) # , marker='→'
        scatterlines!(t, p_xs, label=L"\vec{p}_x", color = Cycled(2)) # , marker='→'
        ax2 = Axis(fig[2, 1], ylabel="Y component")
        # hidexdecorations!(ax2, grid=false)
        scatterlines!(t, r_ys, label=L"\vec{r}_y", color = Cycled(1)) # , marker='↑'
        scatterlines!(t, p_ys, label=L"\vec{p}_y", color = Cycled(2)) # , marker='↑'
        linkxaxes!(ax1, ax2)
        Legend(fig[1, 1], ax1, tellwidth=false, tellheight=false, nbanks=2, valign=:top)
        Legend(fig[2, 1], ax2, tellwidth=false, tellheight=false, nbanks=2, valign=:top)
    # end
    # begin  # Plotting orbits
        # fig = Figure()
        ax = Axis(fig[1:2, 2], title="Orbit", ylabel="Y coordinate", xlabel="X coordinate", yaxisposition=:right, aspect=DataAspect())
        scatterlines!(r_xs, r_ys, label=L"\vec{r}", color = Cycled(1), markersize=t ./ t[end] .* 10)
        scatterlines!(p_xs, p_ys, label=L"\vec{p}", color = Cycled(2), markersize=t ./ t[end] .* 10)
        Legend(fig[:, 2], ax, tellwidth=false, tellheight=false, halign=:left, valign=:top)
    # end
    # r⃗ = positions .- real.(us_analytical) # Residuals
    # errors = [√sum(r⃗[begin:i].^2 .* h) for i in eachindex(r⃗)]
    
    
    fig |> display
end


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

# Solving a kepler problem with euler method
t₀=0
t₁=2*pi*13^(3/2)
N=10000
t=range(t₀,t₁,N)
x=13
y=0
px=0
py=1/sqrt(x) # p=1/sqrt(x) to get circular orbit.
x₀=[x; y; px; py]



xtable=euler(d_dt_kepler,x₀,t)
# plot_orbit(xtable,t)
plot_orbit2(xtable,t)


# let 
    ω=norm([x,y])
    r_anal=[reshape(cos.(ω*t), 1, :); reshape(sin.(ω*t), 1, :)]
    p_anal=ω*[reshape(-sin.(ω*t), 1, :); reshape(cos.(ω*t), 1, :)]
# end


