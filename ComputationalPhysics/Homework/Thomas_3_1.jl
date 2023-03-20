using LinearAlgebra
using Plots

#PROBLEM SET 2

include("Thomas_3_euler.jl")
include("Thomas_3_ref_function.jl")

#Plotly is called to enable interactive methods
plotly()


#1B: Implement Euler method for 2x2 matrix problem
γ = 0
function system_matrix(γ)
    matrix =     [  0 1;
                   -1 -γ]
    return matrix
end
initial_values = [1;0]

prob1 = euler_integration(system_matrix(γ),initial_values,0,30,0.01)
prob1_ref = reference_function(prob1[3,:],γ)

plot(prob1[3,:],prob1[1,:],labels="numeric")

#1C: analytically find the two L.I. solutions to the ODE
#and create reference solutions
plot!(prob1[3,:],real(prob1_ref[:,1]),labels="analytical")
Plots.ylims!(-1.2,1.2)
DHO = plot!(title = "Damped Harmonic Oscillator")
display(DHO)

step_size = [1 0.1 0.01 0.001]
gamma_vals = [0 0.3 2 4]
timesteps = prob1[3,:]

error_plot = plot()
convergence_plot = plot()
convergence = [0 0;]

for n in 1:4
    numerical = euler_integration(system_matrix(gamma_vals[4]),initial_values,0,30,step_size[n])
    analytic = reference_function(numerical[3,:],gamma_vals[4])
    err = @. abs(numerical[1,:] - real(analytic))
    plot!(error_plot,numerical[3,:],err,yaxis=:log)
    scatter!(convergence_plot,(step_size[n],last(err)),yaxis=:log,xaxis=:log)
end

display(convergence_plot)
display(error_plot)

#Problem 2:

#Note that the variables in this case are *vectors*
N = 500
r = zeros(2,N)
p = zeros(2,N)
#Position vector r initial condition
r[:,1] = [1; 0]

#Momentum vector p initial condition
p[:,1] = [0; 0.9]

time_vector = LinRange(0,10,N)
dt = time_vector[2] - time_vector[1]

for n in 2:N
    r[:,n] = r[:,n-1] + dt * p[:,n-1]
    p[:,n] = p[:,n-1] - dt * r[:,n-1] ./ norm(r[:,n-1])^3
    
end

kepler = plot(r[1,:],r[2,:],aspect_ratio=:equal)
display(kepler)



