using LinearAlgebra 


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


using LaTeXStrings
using GLMakie; Makie.inline!(true)
function plot_orbit2(xtable,t,r_analytic) # Credit Dennis. Changed to fit my implementation -Andi
    
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
        scatterlines!(t, r_xs, label=L"\vec{r}_x", color = Cycled(1), markersize=2) # , marker='→'
        scatterlines!(t, p_xs, label=L"\vec{p}_x", color = Cycled(2), markersize=2) # , marker='→'
        scatterlines!(t, r_analytic[1,:], label=L"\vec{r}_{x,anal}", color = :red, markersize=2) # , marker='→'
        ax2 = Axis(fig[2, 1], ylabel="Y component")
        # hidexdecorations!(ax2, grid=false)
        scatterlines!(t, r_ys, label=L"\vec{r}_y", color = Cycled(1), markersize=2) # , marker='↑'
        scatterlines!(t, p_ys, label=L"\vec{p}_y", color = Cycled(2), markersize=2) # , marker='↑'
        scatterlines!(t, r_analytic[2,:], label=L"\vec{r}_{y,anal}", color = :red, markersize=2) # , marker='→'
        linkxaxes!(ax1, ax2)
        Legend(fig[1, 1], ax1, tellwidth=false, tellheight=false, nbanks=2, valign=:top)
        Legend(fig[2, 1], ax2, tellwidth=false, tellheight=false, nbanks=2, valign=:top)
    # end
    # begin  # Plotting orbits
        # fig = Figure()
        ax = Axis(fig[1:2, 2], title="Orbit", ylabel="Y coordinate", xlabel="X coordinate", yaxisposition=:right, aspect=DataAspect())
        scatterlines!(r_xs, r_ys, label=L"\vec{r}", color = Cycled(1), markersize=t ./ t[end] .* 10)
        scatterlines!(p_xs, p_ys, label=L"\vec{p}", color = Cycled(2), markersize=t ./ t[end] .* 10)
        scatterlines!(r_analytic[1,:], r_analytic[2,:], label=L"\vec{r}_{anal}", color = :red, markersize=t ./ t[end] .* 5)
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

x=3
y=0
px=0
py=1/sqrt(x) # p=1/sqrt(x) to get circular orbit.
x₀=[x; y; px; py]

t₀=0
t₁=(2*pi*x^(3/2))*10 #Normalize to whole round trips
N=10000
t=range(t₀,t₁,N)


xtable=euler(d_dt_kepler,x₀,t)

ω=1/norm([x₀[1];x₀[2]])^(3/2)
r_anal=norm([x₀[1];x₀[2]])*[reshape(cos.(ω*t), 1, :); reshape(sin.(ω*t), 1, :)]
p_anal=norm([x₀[1];x₀[2]])*ω*[reshape(-sin.(ω*t), 1, :); reshape(cos.(ω*t), 1, :)]


residual=norm.(eachcol(xtable[1:2,:]))-norm.(eachcol(r_anal[1:2,:]))
# plot_orbit(xtable,t)
plot_orbit2(xtable,t,r_anal)


# let 


using Plots
plotly() # Interactive plots

# Plot error. NOTE +1E-16, to get ticks on plot. Due to points at zero, so log(0)=-∞
Plots.plot(t,residual, labels="Residual")
Plots.plot!(title = "Residual vs. time, x(t=0)=$x, N=$N")
Plots.plot!(legend=:right)
Plots.xlabel!("t [s]")
error_plot1=Plots.ylabel!("Residual [ ]")

Plots.plot(t,residual, labels="Residual", yaxis=:log)
Plots.ylims!(1E-2,maximum(residual))
# plot!(xticks=(1:10, 1:10), grid=true)
Plots.plot!(title = "Residual vs. time, x(t=0)=$x, N=$N")
Plots.plot!(legend=:right)
Plots.xlabel!("t [s]")
error_plot2=Plots.ylabel!("Residual [ ]")


#  plot(t,residual, labels="Residual", xaxis=:log, yaxis=:log)
#  Plots.ylims!(1E-16,maximum(residual))
#  Plots.xlims!(1E-2,maximum(t))
#  # plot!(xticks=(1:10, 1:10), grid=true)
#  plot!(title = "Residual vs. time, γ=$γ, N=$N")
#  plot!(legend=:right)
#  xlabel!("t [s]")
#  error_plot3=ylabel!("Residual [ ]")

error_plot=Plots.plot(error_plot1, error_plot2, layout = (2,1))
display(error_plot)

h=step(t)
function error_vs_time(residual,h,N)
    err_N=√(sum(h.*(residual[1:N]).^2))
    return err_N
end




# end


