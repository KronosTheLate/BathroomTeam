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

N=5000
t₁=100
t₀=0
t=range(t₀,t₁,N)
h=step(t)


γ=0.15
A=-1*[0 -1;
      1 γ]
x₀=[1;0]

xtable=euler(A,x₀,t)

# Analytic solution
if γ==2
    γ+=2*eps()
end

α=x₀[1]
β=x₀[2]
Φ=√((γ/2)^2-1+0im)
@show A=α/4*γ/Φ+α/2+β
@show B=α-A

s₁=-γ/2+√((γ/2)^2-1+0im)
s₂=-γ/2-√((γ/2)^2-1+0im)
anal=@.A*exp(s₁*t)+B*exp(s₂*t)


# Just plotting. Should be changed for different use case
plot(t,xtable[1,:], labels="Position")
plot!(t,real.(anal),labels="Analytic position")
plot!(t,xtable[2,:], labels="Velocity")
plot!(title = "Damped Harmonic Oscillator")
# Plots.plot!(legend=:topright)
xlabel!("t [s]")
harm_osc=ylabel!("Amplitude [m, m/s]")
display(harm_osc)



##



