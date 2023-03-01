using GLMakie; Makie.inline!(true)
using LinearAlgebra  # for `norm` of vectors

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
    timespan = 10
    N = 10^3
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
    Ns = [round(Int64, 10 * 2^i) for i in 0:16]|>reverse
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
    power = 1; power_prop(h) = h^power
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