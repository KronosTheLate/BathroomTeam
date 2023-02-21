using GLMakie; Makie.inline!(true)
using LinearAlgebra  # for `norm` of vectors
#=
function integrate_euler(state_to_rateofchange_matrix, t₀, t₁, steps, u⃗₀)
    timesteps = range(t₀, t₁, steps)
    h = step(timesteps)
    u⃗s = fill(69.0, length(u⃗₀), length(timesteps))  # For outputting
    u⃗s[:, 1] = u⃗₀  # init
    for i in 2:lastindex(timesteps)
        rate_of_change = state_to_rateofchange_matrix * u⃗s[:, i-1]
        Δu⃗ = rate_of_change * h
        u⃗s[:, i] = u⃗s[:, i-1] + Δu⃗
    end
    return (us = u⃗s, timesteps=timesteps)
end
integrate_euler(args) = integrate_euler(args...)
=#

function integrate_euler(∇::Function, t₀, t₁, steps, u⃗₀)
    timesteps = range(t₀, t₁, steps)
    h = step(timesteps)
    u⃗s = Vector{typeof(u⃗₀)}(undef, length(timesteps))  #fill(69.0, length(u⃗₀), length(timesteps))  # For outputting
    u⃗s[1] = u⃗₀  # init
    for i in 2:lastindex(timesteps)
        rate_of_change = ∇(u⃗s[i-1])
        Δu⃗ = rate_of_change * h
        u⃗s[i] = u⃗s[i-1] .+ Δu⃗
    end
    return (us = u⃗s, timesteps=timesteps)
end
integrate_euler(args) = integrate_euler(args...)


let #? Task 1, plotting single case
    begin
        t0 = 0
        t1 = 10
        N = 500
        γ = 0.15
        if γ == 2
            γ += 2eps()  # Workaround
        end
        A(u⃗) = [0.0 1; -1 -γ] * u⃗
        prob = (A, t0, t1, N, [1.0, 0])
    end

    us, timesteps = integrate_euler(prob)
    h = step(timesteps)
    positions = getindex.(us, 1)
    velocities = getindex.(us, 2)
    fig, ax, plt = lines(timesteps, positions, label="Positions")
    α = prob[end][1]
    β = prob[end][2]
    Φ = √(γ^2/4-1+0im)
    B = α/2 - α/4*γ/Φ - β# (Φ - γ/2)/Φ
    A = 1-B
    
    us_analytical = [A*exp(t*(-γ/2 + √(γ^2/4-1+0im))) + 
                        B*exp(t*(-γ/2 - √(γ^2/4-1+0im))) for t in timesteps
    ]
    r⃗ = positions .- real.(us_analytical) # Residuals
    errors = [√sum(r⃗[begin:i].^2 .* h) for i in eachindex(r⃗)]
    # display(errors)

    lines!(timesteps, us_analytical.|>real, label="Analytical position")
    lines!(timesteps, velocities, label="Velocities")
    axislegend()
    # ax.xscale=log10
    fig |> display
end

##¤ Task 1d
function task_1d(γ, timestep)
    t0 = 0
    t1 = 100
    
    N = round(Int64, (t1-t0)/timestep)
    # println("Desired step: $timestep")
    # println("True step   : $(step(range(0, 100, N)))")
    if γ == 2
        γ += 2eps()  # Workaround
    end
    A(u⃗) = [0.0 1; -1 -γ] * u⃗

    prob = (A, t0, t1, N, [1.0, 0])
    us, timesteps = integrate_euler(prob)
    h = step(timesteps)
    positions = getindex.(us, 1)
    α = prob[end][1]
    β = prob[end][2]
    Φ = √(γ^2/4-1+0im)
    B = α/2 - α/4*γ/Φ - β# (Φ - γ/2)/Φ
    A = 1-B
    us_analytical = [A*exp(t*(-γ/2 + √(γ^2/4-1+0im))) + 
                     B*exp(t*(-γ/2 - √(γ^2/4-1+0im))) for t in timesteps
    ]
    r⃗ = positions .- real.(us_analytical) # Residuals
    errors = [√sum(r⃗[begin:i].^2 .* h) for i in eachindex(r⃗)]

    return (timesteps = timesteps, pos_est=positions, pos_ana=real.(us_analytical), resids=r⃗, errors=errors)
end

with_theme(markersize=5, resolution=(1920÷1.75, 1080÷1.5)) do
    for γ in (0.5, 1, 2, 4) # Un, under, critical, over
        fig = Figure()
        for (i, timestep) in enumerate([1, 1e-1, 1e-2, 1e-3][1:3])
            timesteps, pos_est, pos_ana, resids, errors = task_1d(γ, timestep)
            let #! Linear-scale plot
                ax = Axis(fig[1, i], title="timestep = $timestep")
                scatterlines!(timesteps, pos_est, label="Estimated position")
                scatterlines!(timesteps, resids,  label="Residuals")
                scatterlines!(timesteps, pos_ana, label="Analytical position")
                scatterlines!(timesteps, errors,  label="Error up to timestep")
                Legend(fig[4, 2:3], ax, nbanks=2)
                Label(fig[1, 0], "Lin x & y", rotation=π/2, tellheight=false)
            end
            let #! Log-scale plots
                thresh_y = 1e-20
                thresh_timesteps = 1e-1
                increase_floor(v::AbstractVector, thresh=thresh_y) = replace(x->abs(x) ≥ thresh ? x : thresh, v) 
                logproc(x) = x .|> abs |> increase_floor
                timesteps = increase_floor(timesteps, thresh_timesteps)
                ax = Axis(fig[2, i])

                scatterlines!(timesteps, logproc(pos_est))
                scatterlines!(timesteps, logproc(resids) )
                scatterlines!(timesteps, logproc(pos_ana))
                scatterlines!(timesteps, logproc(errors) )
                # ax.xscale = log10
                ax.yscale = log10
                Label(fig[2, 0], "Lin x, Log y", rotation=π/2, tellheight=false)

                ax = Axis(fig[3, i])
                scatterlines!(timesteps, logproc(pos_est))
                scatterlines!(timesteps, logproc(resids) )
                scatterlines!(timesteps, logproc(pos_ana))
                scatterlines!(timesteps, logproc(errors) )
                ax.yscale = log10
                ax.xscale = log10
                Label(fig[3, 0], "Log x & y", rotation=π/2, tellheight=false)
            end
        end
        Label(fig[4, 1], "γ = $γ", tellwidth=false, textsize=25)
        fig |> display
    end
end

##¤ Task 2
# state_to_rateofchange_matrix goes from const to function of state
# The function is renamed to `∇`, and is applied to the current state

with_theme(resolution = (1920÷2, 1080÷2) ,markersize=5) do
    #? Problem definition
    ∇(u⃗) = [0 1; -1/norm(u⃗)^3 0]*u⃗
    t0 = 0
    t1 = 50
    N = 2000
    p⃗₀ = 0.5
    prob = (∇, t0, t1, N, [[1.0, 0.0], [0, p⃗₀]])

    us, timesteps = integrate_euler(prob)
    h = step(timesteps)
    rs = getindex.(us, 1)
    r_xs = getindex.(rs, 1)
    r_ys = getindex.(rs, 2)

    ps = getindex.(us, 2)
    p_xs = getindex.(ps, 1)
    p_ys = getindex.(ps, 2)

    #? Plotting
    begin  # Plotting components
        fig = Figure()
        ax1 = Axis(fig[1, 1], xlabel="Timestep", ylabel="X component")
        scatterlines!(timesteps, r_xs, label=L"\vec{r}_x", color = Cycled(1)) # , marker='→'
        scatterlines!(timesteps, p_xs, label=L"\vec{p}_x", color = Cycled(2)) # , marker='→'
        ax2 = Axis(fig[2, 1], ylabel="Y component")
        scatterlines!(timesteps, r_ys, label=L"\vec{r}_y", color = Cycled(1)) # , marker='↑'
        scatterlines!(timesteps, p_ys, label=L"\vec{p}_y", color = Cycled(2)) # , marker='↑'
        linkxaxes!(ax1, ax2)
        # Legend(fig[2, 1], ax, tellwidth=false, tellheight=true, nbanks=2)
    end
    begin  # Plotting orbits
        # fig = Figure()
        ax = Axis(fig[1:2, 3], title="Orbit", ylabel="Y coordinate", xlabel="X coordinate", yaxisposition=:right)
        scatterlines!(r_xs, r_ys, label=L"\vec{r}", color = Cycled(1), markersize=timesteps ./ t1 .* 10)
        scatterlines!(p_xs, p_ys, label=L"\vec{p}", color = Cycled(2), markersize=timesteps ./ t1 .* 10)
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
    for timestep in (1, 1e-1, 1e-2, 1e-3)
        counter += 1
        timesteps, pos_est, pos_ana, resids, errors = task_1d(γ, timestep)
        let #! Linear-scale plots
            fig, ax, plt = lines(timesteps, pos_est, label="Estimated position")
            ax.title = "γ = $γ, h = $timestep"
            lines!(timesteps, pos_ana, label="Analytical position")
            lines!(timesteps, resids, label="Residuals")
            lines!(timesteps, errors, label="Error up to timestep")
            # axislegend()
            Legend(fig[1, 2], ax)
            fig |> display
        end
        let #! Log-scale plots
            thresh_y = 1e-20
            increase_floor(v::AbstractVector, thresh=thresh_y) = replace(x->abs(x) ≥ thresh ? x : thresh, v) 
            thresh_timesteps = 1e-1
            timesteps = increase_floor(timesteps, thresh_timesteps)
            fig = Figure()
            ax = Axis(fig[1, 1])
            mytitle = "γ = $γ, h = $timestep"
            # mytitle *= "\nthresh_y = $thresh_y, thresh_x = $thresh_timesteps"
            ax.title = mytitle
            scatterlines!(timesteps, abs.(pos_est) |> increase_floor, label="Estimated position")
            scatterlines!(timesteps, increase_floor(abs.(pos_ana)) , label="Analytical position")
            scatterlines!(timesteps, increase_floor(abs.(resids))  , label="Residuals")
            scatterlines!(timesteps, increase_floor(abs.(errors))  , label="Error up to timestep")
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