using LinearAlgebra # Maybe, to get norm
using Plots
plotly() # Interactive plots
# gr()

# using Symbolics
# @variables y;


A=-1*[0 -1;
      1 2]
display(eigen(A))
# Test for critically damped case. Just to check analytical result.
A=-1*[0 -1;
      1 0]
display(eigen(A))
# Test for undamped case. Just to check analytical result.

##

# Euler function for diff.eq. being a simple matrix
function euler(A,x₀,t)
    x=x₀
    h=step(t)
    xtot=zeros(size(x,1),size(t,1))
    xtot[:,1]=x
    for i in 2:size(t,1)
        dx=A*x
        x+=dx*h
        xtot[:,i]=x
    end
    
    return xtot
end
euler(args)=euler(args...)


# Calculating analytic result for specific damping/'factor γ. See OneNote notes for derivation.
function anal_harm_osc(γ,x₀,t)
    # Analytic solution
    if γ==2 # Workaround to avoid problem at critical damped case. Think it is due to very different to linearly independent solutions.
        γ+=2*eps()
    end

    α=x₀[1]
    β=x₀[2]
    Φ=√((γ/2)^2-1+0im)
    A=α/4*γ/Φ+α/2+β
    B=α-A
    # @show A
    # @show B

    s₁=-γ/2+√((γ/2)^2-1+0im)
    s₂=-γ/2-√((γ/2)^2-1+0im)
    anal=@.A*exp(s₁*t)+B*exp(s₂*t)
    return anal
end 


#Plotting the Harmonic oscillator. Both the simulated and analytic position as well as velocity.
#Furthermore the error is plotted in three different axis scales.
function plot_harm_osc(xtable2,anal,t,γ)
    N=size(t,1)
    # Just plotting data and analytical result
    plot(t,xtable2[1,:], labels="Position")
    plot!(t,real.(anal),labels="Analytic position")
    plot!(t,xtable2[2,:], labels="Velocity")
    plot!(title = "Damped Harmonic Oscillator, γ=$γ, N=$N")
    # plot!(legend=:topright)
    xlabel!("t [s]")
    harm_osc=ylabel!("Amplitude [m, m/s]")
    display(harm_osc)

    # Residual and error calculations
    residual=xtable2[1,:]-real.(anal)

    function error_vs_time(residual,h,N)
        err_N=√(sum(h.*(residual[1:N]).^2))
        return err_N
    end

    E=zeros(size(t,1))
    h=step(t)
    for i in 1:size(t,1)
        # @show typeof(E)
        E[i]=error_vs_time(residual,h,i)
    end
    # display(E)

    # Plot error. NOTE +1E-16, to get ticks on plot. Should not be there, but log(0)=-∞
    plot(t.+1E-16,E.+1E-16, labels="Error")
    plot!(title = "Error vs. time, γ=$γ, N=$N")
    Plots.plot!(legend=:right)
    xlabel!("t [s]")
    error_plot1=ylabel!("Error [ ]")


    plot(t.+1E-16,E.+1E-16, labels="Error", yaxis=:log)
    # plot!(xticks=(1:10, 1:10), grid=true)
    plot!(title = "Error vs. time, γ=$γ, N=$N")
    Plots.plot!(legend=:right)
    xlabel!("t [s]")
    error_plot2=ylabel!("Error [ ]")


    plot(t.+1E-16,E.+1E-16, labels="Error", xaxis=:log, yaxis=:log)
    # plot!(xticks=(1:10, 1:10), grid=true)
    plot!(title = "Error vs. time, γ=$γ, N=$N")
    Plots.plot!(legend=:right)
    xlabel!("t [s]")
    error_plot3=ylabel!("Error [ ]")

    error_plot=plot(error_plot1, error_plot2, error_plot3, layout = (3,1))
    display(error_plot)

end

#Just a concice way to formulate this damped harmonic oscillator problem, and make later code more readable.
function harm_osc_prob(t₀,t₁,N,γ,x₀)
    t=range(t₀,t₁,N)
    A=-1*[0 -1;
        1 γ]
    return (A, x₀, t)
end

## Using Euler method to calculate movement of damped harmonic oscillation

let
    N=1000
    t₁=100
    t₀=0
    t=range(t₀,t₁,N)
    # h=step(t)

    γ=0.6
    A=-1*[0 -1;
        1 γ]
    x₀=[1;0]


    xtable=euler(A,x₀,t)
    anal=anal_harm_osc(γ,x₀,t)
    plot_harm_osc(xtable,anal,t,γ)
end


##

# Test of total method for task 1, also using "harm_osc_prob()" for readablity
let
    t₀=0
    t₁=100
    N=500
    γ=0.6
    x₀=[1;0]
    problem1=harm_osc_prob(t₀,t₁,N,γ,x₀)
    @show problem1
    xtable=euler(problem1)
    anal=anal_harm_osc(γ,x₀,problem1[3])
    plot_harm_osc(xtable,anal,problem1[3],γ)
end

##

# Loop over all the required parameters and plot it.

# This takes surprisingly long time to run.... Some error or need optimization?
# I think the error calculation takes much of the time. Can choose to not calculate error for each point to speed up code run time.
for γ in [0, 0.4, 2, 3]
    for Δt in [1, 0.1, 0.01, 0.001]
        t₀=0; t₁=100; x₀=[1;0]
        N=round(Int, (t₁-t₀)/Δt)
        problem1=harm_osc_prob(t₀,t₁,N,γ,x₀)
        xtable=euler(problem1)
        anal=anal_harm_osc(γ,x₀,problem1[3])
        plot_harm_osc(xtable,anal,problem1[3],γ)
    end
end





## ------------------
## Task 2, BUT firstly making more general Euler method
## ------------------

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

## Testing more general Euler method on task 1

t₀=0
t₁=100
N=1000
t=range(t₀,t₁,N)
γ=0.6
x₀=[1;0]
problem1=harm_osc_prob(t₀,t₁,N,γ,x₀)
Amult(xval) = -1*[0 -1;1 γ] * xval # Credit to Dennis for nice and easy way. No need for functors

xtable=euler(Amult,x₀,t)
anal=anal_harm_osc(γ,x₀,t)
plot_harm_osc(xtable,anal,t,γ)

##

function plot_orbit(xtable2,t)
    # Just plotting data and analytical result
    plot(t,xtable2[1,:], labels="x")
    plot!(t,xtable2[2,:], labels="y")
    plot!(t,xtable2[3,:], labels="px")
    plot!(t,xtable2[4,:], labels="py")
    py0=xtable2[4,1]
    plot!(title = "Orbiting object variables, py_0=$py0")
    # plot!(legend=:topright)
    xlabel!("t [s]")
    orbit_plot=ylabel!("Normalized distance or momentum [ ]")
    display(orbit_plot)

end

## OBS! Be aware of 2 different plotting function, which uses different plotting library (plotly and makie).
using LaTeXStrings
using GLMakie; Makie.inline!(true)
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
        hidexdecorations!(ax2, grid=false)
        scatterlines!(t, r_ys, label=L"\vec{r}_y", color = Cycled(1)) # , marker='↑'
        scatterlines!(t, p_ys, label=L"\vec{p}_y", color = Cycled(2)) # , marker='↑'
        linkxaxes!(ax1, ax2)
        # Legend(fig[2, 1], ax, tellwidth=false, tellheight=true, nbanks=2)
    # end
    # begin  # Plotting orbits
        # fig = Figure()
        ax = Axis(fig[1:2, 3], title="Orbit", ylabel="Y coordinate", xlabel="X coordinate", yaxisposition=:right)
        scatterlines!(r_xs, r_ys, label=L"\vec{r}", color = Cycled(1), markersize=t ./ t[end] .* 10)
        scatterlines!(p_xs, p_ys, label=L"\vec{p}", color = Cycled(2), markersize=t ./ t[end] .* 10)
        Legend(fig[:, 2], ax)#, tellwidth=false, tellheight=true)
    # end
    # r⃗ = positions .- real.(us_analytical) # Residuals
    # errors = [√sum(r⃗[begin:i].^2 .* h) for i in eachindex(r⃗)]
    
    
    fig |> display
end


## Task 2. Implementing Kepler, and plotting with above functions

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

# a=d_dt_kepler([1; 0; 0;1]) # Small test of diff.eq.
# display(a)


# Solving a kepler problem with euler method
t₀=0
t₁=60
N=1000
t=range(t₀,t₁,N)
x=1
y=0
px=0
py=1
x₀=[x; y; px; py]
xtable=euler(d_dt_kepler,x₀,t)
# plot_orbit(xtable,t)
plot_orbit2(xtable,t)


## Solving Kepler problem for all different parameters listed in problem sheet

x=1
y=0
px=0
t₀=0
t₁=60

for py in [0, 0.3, 1, 2]
    x₀=[x; y; px; py]
    for Δt in [1, 0.1, 0.01, 0.001]
        N=round(Int, (t₁-t₀)/Δt)
        t=range(t₀,t₁,N)

        xtable=euler(d_dt_kepler,x₀,t)
        # plot_orbit(xtable,t)
        plot_orbit2(xtable,t)
    end
end



