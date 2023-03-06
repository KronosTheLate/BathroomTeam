using GLMakie; Makie.inline!(true)
using LinearAlgebra  # for `norm` of vectors

"""
    step_euler!(u⃗, ∂ₜ, stepsize)

Given a state u⃗, a function ∂ₜ that gives the 
slope of each state when applied to u⃗, and a 
stepsize, update the state by the euler method.
"""
function step_euler!(∂ₜ, u⃗, stepsize)
    # slopes = ∂ₜ(u⃗)
    # Δu⃗ = stepsize * slopes
    # @. u⃗ = u⃗ + Δu⃗
    u⃗ .+= stepsize * ∂ₜ(u⃗)     # The 3 lines above in one
    return nothing
end

"""
    step_leapfrog!(u⃗, ∂ₜ, stepsize)

Given a state u⃗, a function ∂ₜ that gives the 
slope of each state when applied to u⃗, and a 
stepsize, update the state by the leapfrog method.
This means updating the first half of the states, 
recomputing the slopes, and then updating the rest of the states.
"""
function step_leapfrog!(∂ₜ, u⃗, stepsize)  # Taking some care to remain general for 3D
    u⃗[begin:end÷2] .+= stepsize * ∂ₜ(u⃗)[begin:end÷2]  # Updating position only
    u⃗[end÷2+1:end] .+= stepsize * ∂ₜ(u⃗)[end÷2+1:end]  # u⃗ is partially updated, and used to calculate new slope
    return nothing
end

function step_series(∂ₜ, tspan, N_points, u⃗₀, stepper!)
    output = Matrix{Float64}(undef, length(u⃗₀), N_points)
    u⃗ = u⃗₀  # State. Chanes each step
    output[:, 1] = u⃗

    # tpoints are needed for plotting after, 
    # and to get stepsize
    tpoints = range(0, tspan, N_points)
    stepsize = step(tpoints)
    for i in 2:N_points
        stepper!(∂ₜ, u⃗, stepsize)
        output[:, i] = u⃗
    end
    return (us = output, ts=tpoints)
end

let #? Plotting orbit and error
    #* Problem definition
    function ∂ₜ(u⃗)
        rx, ry, px, py = u⃗
        length_r_cubed = hypot(rx, ry)^3
        return [px, py, -rx/length_r_cubed, -ry/length_r_cubed]  # Return SVector for performance?
    end

    tspan = 2π * 10
    N_points = 10^4
    stepper! = step_euler!
    stepper! = step_leapfrog!
    r⃗₀ = [1.0, 0.0]
    p⃗₀ = [0, 1/√r⃗₀[1]]
    u⃗₀ = vcat(r⃗₀, p⃗₀)

    u⃗_ref(t, Ω=1) = [cos(Ω*t), sin(Ω*t), -sin(Ω*t), cos(Ω*t)]
    #? Getting the solution
    solver_name = typeof(stepper!) == typeof(step_euler!) ? "Euler" : 
                  typeof(stepper!) == typeof(step_leapfrog!) ? "Leapfrog" : 
                  "Solver not found"
    us, ts = step_series(∂ₜ, tspan, N_points, u⃗₀, stepper!)
    h = step(ts)

    rs = @view us[1:2, :]
    ps = @view us[3:4, :]

    #? Analytical reference solution:
    us_ref = hcat([u⃗_ref.(t) for t in ts]...)  #! wring shape
    rs_ref = @view us_ref[1:2, :]
    ps_ref = @view us_ref[3:4, :]
        
    #? Plotting orbits
    fig = Figure(resolution=(1080÷2, 1080÷2))
    ax = Axis(fig[1, 1], aspect=1)
    @views scatterlines!(rs[1, :], rs[2, :], label="$solver_name sol",              markersize=range(2, 10, N_points))
    @views scatterlines!(rs_ref[1, :], rs_ref[2, :], label="Analytical sol", markersize=range(2, 10, N_points))
    axislegend()
    fig |> display
    
    #? Plotting error
    resids_r = rs .- rs_ref

    energies_sim = [1/2 * norm(p)^2 - 1/norm(r) for (r, p) in zip(eachcol(rs), eachcol(ps))]
    energies_ref = [1/2 * norm(p)^2 - 1/norm(r) for (r, p) in zip(eachcol(rs_ref), eachcol(ps_ref))]
    resids_E = abs.(energies_sim .- energies_ref)
    fig = Figure(resolution=(1080÷2, 1080÷2));
    ax = Axis(fig[1, 1]);
    
    lines!(ts,     norm.(eachcol(rs)), label="$solver_name |r⃗|")
    lines!(ts,     norm.(eachcol(rs_ref)), label="Ref |r⃗|")
    lines!(ts, norm.(eachcol(resids_r)), label="||r⃗_ref - r⃗||")
    lines!(ts,      resids_E, label="|E_ref - E|")
    # @show extrema(energies[begin+ind_offset:end] .- minimum(energies))
    # xlims!(1, tspan)
    # ylims!(10^-8, 1.3)
    # ax.yscale=log10
    axislegend(position=(1, 0))
    fig|>display
end




##! Old, before mutating stepper functions
function integrate_euler(∇::Function, timespan, steps, u⃗₀)
    timepoints = range(0, timespan, steps)
    h = step(timepoints)
    u⃗s = Vector{typeof(u⃗₀)}(undef, length(timepoints))
    u⃗s[1] = u⃗₀  # init
    for i in 2:lastindex(timepoints)
        rate_of_change = ∇(u⃗s[i-1])
        Δu⃗ = rate_of_change * h
        u⃗s[i] = u⃗s[i-1] .+ Δu⃗
    end
    return (us = u⃗s, timepoints=timepoints)
end
integrate_euler(args) = integrate_euler(args...)


##¤ d) Leapfrog
function integrate_leapfrog(∇::Function, timespan, steps, u⃗₀)
    timepoints = range(0, timespan, steps)
    h = step(timepoints)
    u⃗s = Vector{typeof(u⃗₀)}(undef, length(timepoints))
    u⃗s[1] = u⃗₀  # init
    for i in 2:lastindex(timepoints)
        franken_states = copy(u⃗s[i-1])  # Init as prev state
        for i in eachindex(franken_states)          # Looping over state variables, updating one at a time
            rate_of_change = ∇(franken_states)[i]
            franken_states[i] += rate_of_change * h
        end
        u⃗s[i] = franken_states
    end
    return (us = u⃗s, timepoints=timepoints)
end
integrate_leapfrog(args) = integrate_leapfrog(args...)

##¤ Task 1 b)
#¤  Run the program with this initial condition up to various end 
#¤  times and with various time step sizes that seem plausible. 

let #? Plotting orbit and error
    #? Problem definition
    ∇(u⃗) = [0 1; -1/norm(u⃗[1])^3 0]*u⃗
    timespan = 2π * 10
    N = 10^4
    # solver = integrate_euler
    solver = integrate_leapfrog

    r⃗_ref(t) = let Ω=1; [cos(Ω*t), sin(Ω*t)];end
    p⃗_ref(t) = let Ω=1; [-sin(Ω*t), cos(Ω*t)];end
    r⃗₀ = [1.0, 0.0]
    p⃗₀ = [0, 1/√r⃗₀[1]]
    u⃗₀ = [r⃗₀, p⃗₀]
    prob = (∇, timespan, N, u⃗₀)

    #? Getting the solution
    solver_name = typeof(solver) == typeof(integrate_euler) ? "Euler" : 
                  typeof(solver) == typeof(integrate_leapfrog) ? "Leapfrog" : "Solver not found"
    us, timepoints = solver(prob)
    h = step(timepoints)
    rs = getindex.(us, 1)
    r_xs = getindex.(rs, 1)
    r_ys = getindex.(rs, 2)

    ps = getindex.(us, 2)
    p_xs = getindex.(ps, 1)
    p_ys = getindex.(ps, 2)

    #? Analytical reference solution:
    rs_ref = r⃗_ref.(timepoints)
    r_xs_ref = getindex.(rs_ref, 1)
    r_ys_ref = getindex.(rs_ref, 2)

    ps_ref = p⃗_ref.(timepoints)

    #? Plotting orbits
    fig = Figure(resolution=(1080÷2, 1080÷2))
    ax = Axis(fig[1, 1], aspect=1)
    scatterlines!(r_xs, r_ys, label="$solver_name sol",              markersize=range(2, 10, N))
    scatterlines!(r_xs_ref, r_ys_ref, label="Analytical sol", markersize=range(2, 10, N))
    axislegend()
    fig |> display

    #? Plotting error
    resids_r = [r_ref .- r for (r, r_ref) in zip(rs, rs_ref)]
    energies_sim = [1/2 * norm(p)^2 - 1/norm(r) for (r, p) in zip(rs, ps)]
    energies_ref = [1/2 * norm(p)^2 - 1/norm(r) for (r, p) in zip(rs_ref, ps_ref)]
    resids_E = abs.(energies_sim .- energies_ref)
    fig = Figure(resolution=(1080÷2, 1080÷2));
    ax = Axis(fig[1, 1]);

    lines!(timepoints,     norm.(rs), label="$solver_name |r⃗|")
    lines!(timepoints,     norm.(rs_ref), label="Ref |r⃗|")
    lines!(timepoints, norm.(resids_r), label="||r⃗_ref - r⃗||")
    lines!(timepoints,      resids_E, label="|E_ref - E|")
    # @show extrema(energies[begin+ind_offset:end] .- minimum(energies))
    xlims!(1, timespan)
    ylims!(10^-8, 1.3)
    ax.yscale=log10
    axislegend(position=(1, 0))
    fig|>display
end

#¤      Plot the error E(t) = ||r(t) − r̃(t)|| on a semilog plot versus 
#¤      t and on a doublelog plot versus h. Is the numerical solution 
#¤      stable? If so, what is the convergence order? 
##¤ c) Repeat the same with the error E 0 (t) = |E(t) − Ẽ(t)| of the total energy E(t) = 1/2 * |p|² − |r|⁻¹ .

function get_error_and_stepsize(problem, solver, error_asfun_of_sol)
    sol = solver(problem)
    err = error_asfun_of_sol(sol)
    return (err=err, h=step(sol.timepoints))
end

@time let  #? Convergence plot for 1b)
    Ns = [round(Int64, 10 * 2^i) for i in 0:15]|>reverse
    # stepsizes = [2.0^i for i in -30:2:-15]
    hs = zeros(length(Ns))
    errs = zeros(length(Ns))
    timespan = 15
    # solver = integrate_euler
    solver = integrate_leapfrog
    solver_name = typeof(solver) == typeof(integrate_euler) ? "Euler" : 
                  typeof(solver) == typeof(integrate_leapfrog) ? "Leapfrog" : "Solver not found"
    
    r⃗₀ = [1.0, 0.0]
    p⃗₀ = [0, 1/√r⃗₀[1]]
    u⃗₀ = [r⃗₀, p⃗₀]
    ∇(u⃗) = [0 1; -1/norm(u⃗[1])^3 0]*u⃗
    r⃗_ref(t) = let Ω=1; [cos(Ω*t), sin(Ω*t)];end
    function error_asfun_of_sol_rms(sol)
        L=2
        rs = getindex.(sol.us, 1)
        rs_analytical = r⃗_ref.(sol.timepoints)
        resids_rs = [norm(r .- r_analytical) for (r, r_analytical) in zip(rs, rs_analytical)]
        return sum(resids_rs .^ L .* step(sol.timepoints))^1/L
    end
    error_asfun_of_sol = error_asfun_of_sol_rms
    for (i, N) in enumerate(Ns)
        prob = (∇, timespan, N, u⃗₀)
        errs[i], hs[i] = get_error_and_stepsize(prob, solver, error_asfun_of_sol)
    end
    power = 2; power_prop(h) = h^power
    fig, ax, _ = lines(hs[begin:end÷2], power_prop.(hs[begin:end÷2]) .* errs[1]/power_prop(hs[1]), 
    label="∝ h^$power", linewidth=6
    )
    scatterlines!(hs, errs, color=Cycled(2))
    ax.xscale=log10
    ax.yscale=log10
    ax.xlabel = "Stepsize"
    ax.ylabel = "Error"
    ax.title = solver_name
    axislegend(position=(0, 1))
    fig|>display
end


##