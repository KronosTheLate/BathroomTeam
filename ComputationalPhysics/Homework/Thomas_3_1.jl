using LinearAlgebra
using Plots

#PROBLEM SET 2

include("Thomas_3_euler.jl")
include("Thomas_3_ref_function.jl")

#Plotly is called to enable interactive methods
plotly()


#1B: Implement Euler method for 2x2 matrix problem
γ = 0.3

matrix =     [  0 1;
                -1 -γ]

initial_values = [1;0]

prob1 = euler_integration(matrix,initial_values,0,100,0.1)
prob1_ref = reference_function(prob1[3,:],γ)

plot(prob1[3,:],prob1[1,:],labels="numeric")

#1C: analytically find the two L.I. solutions to the ODE
#and create reference solutions
plot!(prob1[3,:],real(prob1_ref[:,1]),labels="analytical")
Plots.ylims!(-1.2,1.2)
plot!(title = "Damped Harmonic Oscillator")


err =@. abs(prob1[1,:] - real(prob1_ref))
plot(time,err)
