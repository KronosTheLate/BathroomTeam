using LinearAlgebra # Just to get identity matrix
using Plots
plotly() # Interactive plots

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

N=400
t₁=10
t₀=0
t=range(t₀,t₁,N)
h=step(t)


γ=0.5
A=-1*[0 -1;
      1 γ]
x₀=[1;0]



function euler(A,x₀,t)
    x=x₀
    h=step(t)
    xtot=zeros(size(x,1),size(t,1))
    xtot[:,1]=x
    for i in 1:size(t,1)-1
        dx=A*x
        x+=dx*h
        xtot[:,i+1]=x
    end

    # Just plotting
    Plots.plot(t,xtable[1,:], labels="Position")
    Plots.plot!(t,xtable[2,:], labels="Velocity")
    Plots.plot!(title = "Damped Harmonic Oscillator")
    # Plots.plot!(legend=:topright)
    Plots.xlabel!("t [s]")
    harm_osc=Plots.ylabel!("Amplitude [m, m/s]")
    display(harm_osc)
    
    return xtot
end
euler(args)=euler(args...)


xtable=euler(A,x₀,t)






##



