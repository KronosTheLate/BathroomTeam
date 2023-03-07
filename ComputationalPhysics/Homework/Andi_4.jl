using LinearAlgebra 


# More general Euler method. Used for task 2
function euler(f,x₀,t)
    x=copy(x₀)
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

function RK4(f,x₀,t)
    x=copy(x₀)
    h=step(t)
    xtot=zeros(size(x,1),size(t,1))
    xtot[:,1]=x
    for i in 2:size(t,1)
        k1=f(x)
        k2=f(x+h*k1/2)
        k3=f(x+h*k2/2)
        k4=f(x+h*k3)
        x+=1/6*(k1+2*k2+2*k3+k4)*h
        xtot[:,i]=x
    end
    
    return xtot
end
RK4(args)=RK4(args...)



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
        SizeForMarker=t ./ t[end] .* 5
        SizeForMarker=collect(SizeForMarker)
        SizeForMarker[end]*=4
        scatterlines!(r_xs, r_ys, label=L"\vec{r}", color = Cycled(1), markersize=SizeForMarker)
        scatterlines!(p_xs, p_ys, label=L"\vec{p}", color = Cycled(2), markersize=SizeForMarker)
        scatterlines!(r_analytic[1,:], r_analytic[2,:], label=L"\vec{r}_{anal}", color = :red, markersize=SizeForMarker./2)
        Legend(fig[:, 2], ax, tellwidth=false, tellheight=false, halign=:left, valign=:top)
    # end
    # r⃗ = positions .- real.(us_analytical) # Residuals
    # errors = [√sum(r⃗[begin:i].^2 .* h) for i in eachindex(r⃗)]
    
    
    fig |> display
end


using Plots
plotly() # Interactive plots
function plot_residuals(t,residual, residual_E,x,N=size(t,1))

    # Plot error. NOTE +1E-16, to get ticks on plot. Due to points at zero, so log(0)=-∞
    Plots.plot(t,residual, labels="Position residual")
    Plots.plot!(title = "Residual vs. time, x(t=0)=$x, N=$N")
    Plots.plot!(legend=:right)
    Plots.xlabel!("t [s]")

    Plots.plot!(t,residual_E, labels="Energy residual")
    residual1=Plots.ylabel!("Residual [ ]")

    Plots.plot(t,residual, labels="Position residual", yaxis=:log)
    Plots.ylims!(1E-8,maximum(residual))
    # plot!(xticks=(1:10, 1:10), grid=true)
    Plots.plot!(title = "Residual vs. time, x(t=0)=$x, N=$N")
    Plots.plot!(legend=:right)
    Plots.xlabel!("t [s]")

    Plots.plot!(t,residual_E, labels="Energy residual")
    residual2=Plots.ylabel!("Residual [ ]")


    #  plot(t,residual, labels="Residual", xaxis=:log, yaxis=:log)
    #  Plots.ylims!(1E-16,maximum(residual))
    #  Plots.xlims!(1E-2,maximum(t))
    #  # plot!(xticks=(1:10, 1:10), grid=true)
    #  plot!(title = "Residual vs. time, γ=$γ, N=$N")
    #  plot!(legend=:right)
    #  xlabel!("t [s]")
    #  residual3=ylabel!("Residual [ ]")

    residual_plot=Plots.plot(residual1, residual2, layout = (2,1))
    display(residual_plot)
end

function plot_error_vs_h(hs, errors, errors_E, t₁, assumed_order, E_assumed_order=assumed_order)
    t1_time=round(t₁)
    Plots.scatter(hs,errors, labels="Position error", markershape=:cross, yaxis=:log, xaxis=:log)
    Plots.plot!(title = "Error vs. stepsize, t₁=$t1_time")
    Plots.plot!(legend=:bottomright)
    Plots.xlabel!("h [ ]")
    error_vs_h=Plots.ylabel!("Error [ ]")
    f(x)=2*errors[1]*(x)^assumed_order
    Plots.plot!(f, label="Theoretical (position) order: $assumed_order")

    # E_assumed_order=assumed_order
    Plots.scatter!(hs,errors_E, labels="Energy error", markershape=:cross)
    h(x)=6*errors_E[1]*(x)^E_assumed_order
    Plots.plot!(h, label="Theoretical (energy) order: $E_assumed_order")

    display(error_vs_h)
end

function error_vs_time(residual,h,N)
    err_N=√(sum(h.*(residual[1:N]).^2))
    return err_N
end

##


# Vectors for error calculations, and convergence estimation.
hs=[]
errors=[]
errors_E=[]

# Initial conditions
x=3
y=0
px=0
py=1/sqrt(x) # p=1/sqrt(x) to get circular orbit.
x₀=[x; y; px; py]


# Time interval
t₀=0
t₁=(2*pi*x^(3/2))*2 #Normalize to whole round trips
N=10000

# for t₁ in (2*pi*x^(3/2))*[0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4]
# Loop over N
for N in [10, 20, 50, 100, 200, 500, 1000, 2000, 4000, 8000]#, 16000, 32000, 64000]
t=range(t₀,t₁,N)

# xtable=euler(d_dt_kepler,x₀,t)

xtable=leap_frog(d_dt_kepler,x₀,t)
# display(x₀)


# xtable=RK4(d_dt_kepler,x₀,t)


ω=1/norm([x₀[1];x₀[2]])^(3/2)
r_anal=norm([x₀[1];x₀[2]])*[reshape(cos.(ω*t), 1, :); reshape(sin.(ω*t), 1, :)]
p_anal=norm([x₀[1];x₀[2]])*ω*[reshape(-sin.(ω*t), 1, :); reshape(cos.(ω*t), 1, :)]


residual=norm.(eachcol(xtable[1:2,:]))-norm.(eachcol(r_anal)) #.-norm([x₀[1];x₀[2]])
    # Note that i have flipped the order of calculated and analytical/true compared to problem set.
# plot_orbit(xtable,t)
plot_orbit2(xtable,t,r_anal)


E_num=norm.(eachcol(xtable[3:4,:])).*0.5-norm.(eachcol(xtable[1:2,:])).^(-1)
E_anal=norm.(eachcol(p_anal)).*0.5-norm.(eachcol(r_anal)).^(-1)
residual_E=E_num-E_anal
    # Note that i have flipped the order of calculated and analytical/true compared to problem set.


plot_residuals(t,residual, residual_E ,x)

h=step(t)

error_at_h=error_vs_time(residual,h,N)
error_E_at_h=error_vs_time(residual_E,h,N)

append!(hs,h)
append!(errors,error_at_h)
append!(errors_E,error_E_at_h)

end

##

assumed_order=1
# assumed_order=2
# assumed_order=4
plot_error_vs_h(hs, errors, errors_E, t₁, assumed_order, assumed_order)




