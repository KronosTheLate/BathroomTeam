using GLMakie; Makie.inline!(true)
using LinearAlgebra  # for `norm` of vectors
#! I replaced timestep with timepoint for all occurences without checking
#=
function integrate_euler(state_to_rateofchange_matrix, t₀, t₁, steps, u⃗₀)
    timepoints = range(t₀, t₁, steps)
    h = step(timepoints)
    u⃗s = fill(69.0, length(u⃗₀), length(timepoints))  # For outputting
    u⃗s[:, 1] = u⃗₀  # init
    for i in 2:lastindex(timepoints)
        rate_of_change = state_to_rateofchange_matrix * u⃗s[:, i-1]
        Δu⃗ = rate_of_change * h
        u⃗s[:, i] = u⃗s[:, i-1] + Δu⃗
    end
    return (us = u⃗s, timepoints=timepoints)
end
integrate_euler(args) = integrate_euler(args...)
=#

function integrate_euler(∇::Function, timespan, steps, u⃗₀)
    timepoints = range(0, timespan, steps)
    h = step(timepoints)
    u⃗s = Vector{typeof(u⃗₀)}(undef, length(timepoints))  #fill(69.0, length(u⃗₀), length(timepoints))  # For outputting
    u⃗s[1] = u⃗₀  # init
    for i in 2:lastindex(timepoints)
        rate_of_change = ∇(u⃗s[i-1])
        Δu⃗ = rate_of_change * h
        u⃗s[i] = u⃗s[i-1] .+ Δu⃗
    end
    return (us = u⃗s, timepoints=timepoints)
end
integrate_euler(args) = integrate_euler(args...) 

let #? Task 1, plotting single case
    begin
        timespan = 50
        N = 100
        γ = 6
        if γ == 2
            γ += 2eps()  #! Dirty, dirty lie
        end
        ∇(u⃗) = [0.0 1; -1 -γ] * u⃗
        prob = (∇, timespan, N, [1.0, 0])
    end

    us, timepoints = integrate_euler(prob)
    h = step(timepoints)
    positions = getindex.(us, 1)
    velocities = getindex.(us, 2)
    fig, ax, plt = lines(timepoints, positions, label="Positions")
    α = prob[end][1]
    β = prob[end][2]
    Φ = √(γ^2/4-1+0im)
    B = α/2 - α/4*γ/Φ - β# (Φ - γ/2)/Φ
    A = 1-B
    
    us_analytical = [A*exp(t*(-γ/2 + √(γ^2/4-1+0im))) + 
                        B*exp(t*(-γ/2 - √(γ^2/4-1+0im))) for t in timepoints
    ]
    r⃗ = positions .- real.(us_analytical) # Residuals
    errors = [√sum(r⃗[begin:i].^2 .* h) for i in eachindex(r⃗)]
    # display(errors)

    lines!(timepoints, us_analytical.|>real, label="Analytical position")
    lines!(timepoints, velocities, label="Velocities")
    axislegend()
    # ax.xscale=log10
    fig |> display
end

##¤ Task 1d
function task_1d(γ, timepoint)
    timespan = 100
    
    N = round(Int64, (timespan)/timepoint)
    # println("Desired step: $timepoint")
    # println("True step   : $(step(range(0, 100, N)))")
    if γ == 2
        γ += 2eps()  # Workaround
    end
    A(u⃗) = [0.0 1; -1 -γ] * u⃗

    prob = (A, timespan, N, [1.0, 0])
    us, timepoints = integrate_euler(prob)
    h = step(timepoints)
    positions = getindex.(us, 1)
    α = prob[end][1]
    β = prob[end][2]
    Φ = √(γ^2/4-1+0im)
    B = α/2 - α/4*γ/Φ - β# (Φ - γ/2)/Φ
    A = 1-B
    us_analytical = [A*exp(t*(-γ/2 + √(γ^2/4-1+0im))) + 
                     B*exp(t*(-γ/2 - √(γ^2/4-1+0im))) for t in timepoints
    ]
    r⃗ = positions .- real.(us_analytical) # Residuals
    errors = [√sum(r⃗[begin:i].^2 .* h) for i in eachindex(r⃗)]

    return (timepoints = timepoints, pos_est=positions, pos_ana=real.(us_analytical), resids=r⃗, errors=errors)
end

with_theme(markersize=5, resolution=(1920÷1.75, 1080÷1.5)) do
    for γ in (0.5, 1, 2, 4) # Un, under, critical, over
        fig = Figure()
        for (i, timepoint) in enumerate([1, 1e-1, 1e-2, 1e-3][1:3])
            timepoints, pos_est, pos_ana, resids, errors = task_1d(γ, timepoint)
            let #! Linear-scale plot
                ax = Axis(fig[1, i], title="timepoint = $timepoint")
                scatterlines!(timepoints, pos_est, label="Estimated position")
                scatterlines!(timepoints, resids,  label="Residuals")
                scatterlines!(timepoints, pos_ana, label="Analytical position")
                scatterlines!(timepoints, errors,  label="Error up to timepoint")
                Legend(fig[4, 2:3], ax, nbanks=2)
                Label(fig[1, 0], "Lin x & y", rotation=π/2, tellheight=false)
            end
            let #! Log-scale plots
                thresh_y = 1e-20
                thresh_timepoints = 1e-1
                increase_floor(v::AbstractVector, thresh=thresh_y) = replace(x->abs(x) ≥ thresh ? x : thresh, v) 
                logproc(x) = x .|> abs |> increase_floor
                timepoints = increase_floor(timepoints, thresh_timepoints)
                ax = Axis(fig[2, i])

                scatterlines!(timepoints, logproc(pos_est))
                scatterlines!(timepoints, logproc(resids) )
                scatterlines!(timepoints, logproc(pos_ana))
                scatterlines!(timepoints, logproc(errors) )
                # ax.xscale = log10
                ax.yscale = log10
                Label(fig[2, 0], "Lin x, Log y", rotation=π/2, tellheight=false)

                ax = Axis(fig[3, i])
                scatterlines!(timepoints, logproc(pos_est))
                scatterlines!(timepoints, logproc(resids) )
                scatterlines!(timepoints, logproc(pos_ana))
                scatterlines!(timepoints, logproc(errors) )
                ax.yscale = log10
                ax.xscale = log10
                Label(fig[3, 0], "Log x & y", rotation=π/2, tellheight=false)
            end
        end
        Label(fig[4, 1], "γ = $γ", tellwidth=false, textsize=25)
        fig |> display
    end
end

let #! Very manual determining of convergence order
    γ = 2#0.15
    timepoints = [1.1 .^ n for n in -45:20]
    errors = similar(timepoints)
    for (i, timepoint) in enumerate(timepoints)
        sol = task_1d(γ, timepoint)
        errors[i] = sol.errors[end]
    end
    fig, ax, _ = scatterlines(timepoints, errors)
    ax.xscale=log10
    ax.yscale=log10
    ax.xlabel = "Stepsize"
    ax.ylabel = "Error"

    # lines!(timepoints[begin:begin+10], x->10x^1, color=Cycled(2))
    lines!(timepoints[begin:end÷2], x->0.5x^1, color=Cycled(2))
    fig|>display
end

##¤ Task 2
# state_to_rateofchange_matrix goes from const to function of state
# The function is renamed to `∇`, and is applied to the current state

with_theme(resolution = (1920÷2, 1080÷2) ,markersize=5) do
    #? Problem definition
    ∇(u⃗) = [0 1; -1/norm(u⃗[1])^3 0]*u⃗
    timespan = 2π
    N = 10^4

    r⃗₀ = [1.0, 0.0]
    p₀_y = 0.0
    p⃗₀ = [0, p₀_y]
    prob = (∇, timespan, N, [r⃗₀, p⃗₀])

    us, timepoints = integrate_euler(prob)
    h = step(timepoints)
    rs = getindex.(us, 1)
    r_xs = getindex.(rs, 1)
    r_ys = getindex.(rs, 2)

    ps = getindex.(us, 2)
    p_xs = getindex.(ps, 1)
    p_ys = getindex.(ps, 2)

    #? Plotting
    begin  # Plotting components
        fig = Figure()
        ax1 = Axis(fig[1, 1], title="r⃗₀ = $r⃗₀\np⃗₀ = $p⃗₀\nN = $N", xlabel="Timepoint", ylabel="X component",)
        scatterlines!(timepoints, r_xs, label=L"\vec{r}_x", color = Cycled(1)) # , marker='→'
        scatterlines!(timepoints, p_xs, label=L"\vec{p}_x", color = Cycled(2)) # , marker='→'
        
        # scatterlines!(timepoints, (r_xs.|>abs) .+ 1e-3, label=L"\vec{r}_x", color = Cycled(1)) # , marker='→'
        # scatterlines!(timepoints, (p_xs.|>abs) .+ 1e-3, label=L"\vec{p}_x", color = Cycled(2)) # , marker='→'
        # ax1.yscale = log10

        ax2 = Axis(fig[2, 1], ylabel="Y component")
        hidexdecorations!(ax2, grid=false)
        scatterlines!(timepoints, r_ys, label=L"\vec{r}_y", color = Cycled(1)) # , marker='↑'
        scatterlines!(timepoints, p_ys, label=L"\vec{p}_y", color = Cycled(2)) # , marker='↑'
        linkxaxes!(ax1, ax2)
        # Legend(fig[2, 1], ax, tellwidth=false, tellheight=true, nbanks=2)
    end
    begin  # Plotting orbits
        # fig = Figure()
        ax = Axis(fig[1:2, 3], title="Orbit", ylabel="Y coordinate", xlabel="X coordinate", yaxisposition=:right)
        scatterlines!(r_xs, r_ys, label=L"\vec{r}", color = Cycled(1), markersize=timepoints ./ timespan .* 10)
        scatterlines!(p_xs, p_ys, label=L"\vec{p}", color = Cycled(2), markersize=timepoints ./ timespan .* 10)
        Legend(fig[:, 2], ax)#, tellwidth=false, tellheight=true)
    end
    # r⃗ = positions .- real.(us_analytical) # Residuals
    # errors = [√sum(r⃗[begin:i].^2 .* h) for i in eachindex(r⃗)]
    
    
    fig |> display
end



#=
with_theme(markersize=5) do
counter = 0
for γ in (0.5, 1, 2, 4) # Un, under, critical, over
    for timepoint in (1, 1e-1, 1e-2, 1e-3)
        counter += 1
        timepoints, pos_est, pos_ana, resids, errors = task_1d(γ, timepoint)
        let #! Linear-scale plots
            fig, ax, plt = lines(timepoints, pos_est, label="Estimated position")
            ax.title = "γ = $γ, h = $timepoint"
            lines!(timepoints, pos_ana, label="Analytical position")
            lines!(timepoints, resids, label="Residuals")
            lines!(timepoints, errors, label="Error up to timepoint")
            # axislegend()
            Legend(fig[1, 2], ax)
            fig |> display
        end
        let #! Log-scale plots
            thresh_y = 1e-20
            increase_floor(v::AbstractVector, thresh=thresh_y) = replace(x->abs(x) ≥ thresh ? x : thresh, v) 
            thresh_timepoints = 1e-1
            timepoints = increase_floor(timepoints, thresh_timepoints)
            fig = Figure()
            ax = Axis(fig[1, 1])
            mytitle = "γ = $γ, h = $timepoint"
            # mytitle *= "\nthresh_y = $thresh_y, thresh_x = $thresh_timepoints"
            ax.title = mytitle
            scatterlines!(timepoints, abs.(pos_est) |> increase_floor, label="Estimated position")
            scatterlines!(timepoints, increase_floor(abs.(pos_ana)) , label="Analytical position")
            scatterlines!(timepoints, increase_floor(abs.(resids))  , label="Residuals")
            scatterlines!(timepoints, increase_floor(abs.(errors))  , label="Error up to timepoint")
            Legend(fig[1, 2], ax)
            ax.yscale = log10
            fig |> display
            ax.xscale = log10
            fig |> display
        end
        @info "Finnished with plot $counter"
    end
end
end
=#