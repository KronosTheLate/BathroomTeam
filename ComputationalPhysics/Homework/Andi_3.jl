using LinearAlgebra # Just to get identity matrix
using Plots
# plotly() # Interactive plots
gr()

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

##

function anal_harm_osc(γ,x₀,t)
    # Analytic solution
    if γ==2
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




function plot_harm_osc(xtable2,anal,t)
    # Just plotting data and analytical result
    plot(t,xtable2[1,:], labels="Position")
    plot!(t,real.(anal),labels="Analytic position")
    plot!(t,xtable2[2,:], labels="Velocity")
    plot!(title = "Damped Harmonic Oscillator")
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

    # Plot error
    plot(t,E, labels="Error")
    plot!(title = "Error vs. time")
    Plots.plot!(legend=:right)
    xlabel!("t [s]")
    error_plot1=ylabel!("Error [ ]")
    # display(error_plot)


    plot(t,E, labels="Error", yaxis=:log)
    # plot!(xticks=(1:10, 1:10), grid=true)
    plot!(title = "Error vs. time")
    Plots.plot!(legend=:right)
    xlabel!("t [s]")
    error_plot2=ylabel!("Error [ ]")
    # display(error_plot)


    plot(t,E, labels="Error", xaxis=:log, yaxis=:log)
    # plot!(xticks=(1:10, 1:10), grid=true)
    plot!(title = "Error vs. time")
    Plots.plot!(legend=:right)
    xlabel!("t [s]")
    error_plot3=ylabel!("Error [ ]")

    error_plot=plot(error_plot1, error_plot2, error_plot3, layout = (3,1))
    display(error_plot)

end

function harm_osc_prob(t₀,t₁,N,γ1,x₀)
    t=range(t₀,t₁,N)
    A=-1*[0 -1;
        1 γ1]
    return (A, x₀, t)
end

##

let
    N=500
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
    plot_harm_osc(xtable,anal,t)
end


##

# Test of total method for task 1
let
    t₀=0
    t₁=100
    N=500
    γ=0.8
    x₀=[1;0]
    problem1=harm_osc_prob(t₀,t₁,N,γ,x₀)
    @show problem1
    xtable=euler(problem1)
    anal=anal_harm_osc(γ,x₀,problem1[3])
    plot_harm_osc(xtable,anal,problem1[3])
end

##

# Loop over all the required parameters and plot it.

# for Δt in [1, 0.1, 0.001, 0.0001]
#     t₀=0; t₁=100; γ=0.7; x₀=[1;0]
#     N=round(Int, (t₁-t₀)/Δt)
#     problem1=harm_osc_prob(t₀,t₁,N,γ,x₀)
#     xtable=euler(problem1)
#     anal=anal_harm_osc(γ,x₀,problem1[3])
#     plot_harm_osc(xtable,anal,problem1[3])
# end


