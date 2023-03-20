using LinearAlgebra
using Plots

#PROBLEM SET 3

include("Thomas_4_euler.jl")
include("Thomas_4_system_func.jl")
include("Thomas_4_leapfrog.jl")
include("Thomas_4_RK4.jl")

#1A: Initial conditions for circular orbit
rp = [1; 0; 0; 1;]

#1A: reference energy of system (the system has energy conservation)
reference_energy = norm(rp[3:4])^2 / 2 - 1 / norm(rp[1:2])

#1B: Euler method

Tmax = 30
step_size = [0.5 0.2 0.1 0.05 0.02 0.01 0.005 0.002 0.001]

orbit_plot          = plot()
error_plot          = plot()
convergence_plot    = plot()

for n in 1:length(step_size)
    #numerical = euler_integration(system_func,rp,0,Tmax,step_size[n])
    #numerical = leapfrog_integration(system_func,rp,0,Tmax,step_size[n])
    numerical = RK4(system_func,rp,0,Tmax,step_size[n])
    plot!(orbit_plot,numerical[1,:],numerical[2,:],aspect_ratio=:equal)
    print(n)

    #energy of system to get error
    energy =@. (abs(numerical[3,:])^2 + abs(numerical[4,:])^2) / 2 - 1/sqrt(abs(numerical[1,:])^2 + abs(numerical[2,:])^2) 
    err = @. abs(energy - reference_energy)
    plot!(error_plot,numerical[5,:],err,yaxis=:log)

    #convergence
    scatter!(convergence_plot,(step_size[n],last(err)),yaxis=:log,xaxis=:log)
end

#plot!(convergence_plot,step_size,step_size./step_size[1],yaxis=:log,xaxis=:log)

display(orbit_plot)
display(error_plot)
display(convergence_plot)




