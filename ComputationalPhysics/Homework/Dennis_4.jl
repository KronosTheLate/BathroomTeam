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

##¤ Task 1
##¤ b)  Run the program with this initial condition up to various end 
#¤      times and with various time step sizes that seem plausible. 

let #? Plotting orbit and error
    #? Problem definition
    ∇(u⃗) = [0 1; -1/norm(u⃗[1])^3 0]*u⃗
    timespan = 100
    N = 10^4

    r⃗_ref(t) = let Ω=1; [cos(Ω*t), sin(Ω*t)];end
    p⃗_ref(t) = let Ω=1; [-sin(Ω*t), cos(Ω*t)];end
    r⃗₀ = [0.50, 0.0]
    p⃗₀ = [0, 1/√r⃗₀[1]]
    u⃗₀ = [r⃗₀, p⃗₀]
    prob = (∇, timespan, N, u⃗₀)

    #? Getting the solution
    us, timepoints = integrate_euler(prob)
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
    scatterlines!(r_xs, r_ys, label="Euler sol",              markersize=range(2, 10, N))
    scatterlines!(r_xs_ref, r_ys_ref, label="Analytical sol", markersize=range(2, 10, N))
    axislegend()
    # fig |> display

    #? Plotting error
    resids_r = [r_ref .- r for (r, r_ref) in zip(rs, rs_ref)]
    energies_sim = [1/2 * norm(p)^2 - 1/norm(r) for (r, p) in zip(rs, ps)]
    energies_ref = [1/2 * norm(p)^2 - 1/norm(r) for (r, p) in zip(rs_ref, ps_ref)]
    resids_E = abs.(energies_sim .- energies_ref)
    # @show extrema(energies)
    fig = Figure(resolution=(1080÷2, 1080÷2));
    ax = Axis(fig[1, 1]);

    lines!(timepoints,     norm.(rs), label="Euler |r⃗|")
    lines!(timepoints,     norm.(rs_ref), label="Ref |r⃗|")
    lines!(timepoints, norm.(resids_r), label="||r⃗_ref - r⃗||")
    lines!(timepoints,      resids_E, label="|E_ref - E|")
    # @show extrema(energies[begin+ind_offset:end] .- minimum(energies))
    # xlims!(1, timespan)
    ylims!(2^-8, maximum(norm.(resids_r)))
    ax.yscale=log2
    axislegend(position=(1, 0))
    # fig|>display
end

#¤      Plot the error E(t) = ||r(t) − r̃(t)|| on a semilog plot versus 
#¤      t and on a doublelog plot versus h. Is the numerical solution 
#¤      stable? If so, what is the convergence order? 
##¤ c) Repeat the same with the error E 0 (t) = |E(t) − Ẽ(t)| of the total energy E(t) = 1/2 * |p|² − |r|⁻¹ .
#! Error is not a function of time, residuals is. I do not get the semilog plot
#! The convergence order is 0.03 - abissmal.

function get_error_and_stepsize(problem, solver, error_asfun_of_sol)
    sol = solver(problem)
    err = error_asfun_of_sol(sol)
    return (err=err, h=step(sol.timepoints))
end

let  #? Convergence plot for 1b)
    stepsizes = [1.1^i for i in -70:-0]
    hs = zeros(length(stepsizes))
    errs = zeros(length(stepsizes))
    timespan = 100
    ∇(u⃗) = [0 1; -1/norm(u⃗[1])^3 0]*u⃗
    r⃗₀ = [2.0, 0.0]
    p⃗₀ = [0, 1/√r⃗₀[1]]
    u⃗₀ = [r⃗₀, p⃗₀]
    r⃗_ref(t) = let Ω=1; [cos(Ω*t), sin(Ω*t)];end
    function error_asfun_of_sol_rms(sol)
        L=10
        rs = getindex.(sol.us, 1)
        rs_analytical = r⃗_ref.(sol.timepoints)
        resids_rs = [norm(r .- r_analytical) for (r, r_analytical) in zip(rs, rs_analytical)]
        return sum(resids_rs .^ L .* step(sol.timepoints))^1/L
    end
    error_asfun_of_sol = error_asfun_of_sol_rms
    for (i, stepsize) in enumerate(stepsizes)
        N = round(Int64, timespan/stepsize)
        prob = (∇, timespan, N, u⃗₀)
        errs[i], hs[i] = get_error_and_stepsize(prob, integrate_euler, error_asfun_of_sol)
    end
    power = 0.03; power_prop(h) = h^power
    fig, ax, _ = lines(hs[begin:end÷2], power_prop.(hs[begin:end÷2]) .* errs[1]/power_prop(hs[1]), 
    label="∝ h^$power", linewidth=6
    )
    scatterlines!(hs, errs, color=Cycled(2))
    ax.xscale=log2
    ax.yscale=log2
    ax.xlabel = "Stepsize"
    ax.ylabel = "Error"
    axislegend(position=(0, 1))
    fig|>display
end

##¤ d) Leapfrog
#! INCOMPLETE
function integrate_leapfrog(∇::Function, timespan, steps, u⃗₀)
    timepoints = range(0, timespan, steps)
    h = step(timepoints)
    u⃗s = Vector{typeof(u⃗₀)}(undef, length(timepoints))
    u⃗s[1] = u⃗₀  # init
    for i in 2:lastindex(timepoints)
        prev_state = u⃗s[i-1]

        rate_of_change = ∇(u⃗s[i-1])
        Δu⃗ = rate_of_change * h
        u⃗s[i] = u⃗s[i-1] .+ Δu⃗
    end
    return (us = u⃗s, timepoints=timepoints)
end
integrate_leapfrog(args) = integrate_leapfrog(args...)
##
1
# begin  # Plotting components
#     ax1 = Axis(fig[1, 1], title="r⃗₀ = $r⃗₀\np⃗₀ = $p⃗₀\nN = $N", xlabel="Timepoint", ylabel="X component",)
#     scatterlines!(timepoints, r_xs, label=L"\vec{r}_x", color = Cycled(1)) # , marker='→'
#     scatterlines!(timepoints, p_xs, label=L"\vec{p}_x", color = Cycled(2)) # , marker='→'
    
#     # scatterlines!(timepoints, (r_xs.|>abs) .+ 1e-3, label=L"\vec{r}_x", color = Cycled(1)) # , marker='→'
#     # scatterlines!(timepoints, (p_xs.|>abs) .+ 1e-3, label=L"\vec{p}_x", color = Cycled(2)) # , marker='→'
#     # ax1.yscale = log10

#     ax2 = Axis(fig[2, 1], ylabel="Y component")
#     hidexdecorations!(ax2, grid=false)
#     scatterlines!(timepoints, r_ys, label=L"\vec{r}_y", color = Cycled(1)) # , marker='↑'
#     scatterlines!(timepoints, p_ys, label=L"\vec{p}_y", color = Cycled(2)) # , marker='↑'
#     linkxaxes!(ax1, ax2)
# end

# begin  # Plotting orbits
#     ax = Axis(fig[1:2, 3], title="Orbit", ylabel="Y coordinate", xlabel="X coordinate", yaxisposition=:right)
#     scatterlines!(r_xs, r_ys, label=L"\vec{r}", color = Cycled(1), markersize=timepoints ./ timespan .* 10)
#     scatterlines!(p_xs, p_ys, label=L"\vec{p}", color = Cycled(2), markersize=timepoints ./ timespan .* 10)
#     Legend(fig[:, 2], ax)#, tellwidth=false, tellheight=true)
# end
