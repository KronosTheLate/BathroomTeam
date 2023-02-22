using LinearAlgebra
using Plots

include("Thomas_3_euler.jl")

#Plotly is called to enable interactive methods
plotly()

matrix = [  0 -1;
            1  3]

initial_values = [1;2]

prob1 = euler_integration(matrix,initial_values,0,10,0.1)

plot(prob1[3,:],prob1[1,:])

##

