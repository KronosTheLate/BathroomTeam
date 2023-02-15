using GLMakie; Makie.inline!(true)
function integrate(t₀, t₁, steps, state_to_rateofchange_matrix, u⃗₀)
    timesteps = range(t₀, t₁, steps)
    h = step(timesteps)
    u⃗s = fill(69.0, length(u⃗₀), length(timesteps))  # For outputting
    u⃗s[:, 1] = u⃗₀  # init
    for i in 2:lastindex(timesteps)
        Δu⃗ = state_to_rateofchange_matrix * u⃗s[:, i-1] * h
        u⃗s[:, i] = u⃗s[:, i-1] + Δu⃗
    end
    return (us = u⃗s, timesteps=timesteps, h=h)
end
integrate(args) = integrate(args...)

##
t0 = 0
t1 = 100
# desired_hs = (1, 0.1, 0.01, 0.001)

N = 5000
γ = 0.15
if γ == 2
    γ += 2eps()  # Workaround
end
A = [0 1; -1 -γ]

prob1 = (t0, t1, N, A, [1, 0])

# function euler(prob)
let prob = prob1
    us, timesteps, h = integrate(prob)
    step_ = step(timesteps)
    # us|>display
    positions = us[1, :]
    velocities = us[2, :]
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
# euler(prob1)

function task_2d(γ, timestep)
    t0 = 0
    t1 = 100
    
    N = round(Int64, (t1-t0)/timestep)
    println("Desired step: $timestep")
    println("True step   : $(step(range(0, 100, N)))")
    if γ == 2
        γ += 2eps()  # Workaround
    end
    A = [0 1; -1 -γ]

    prob = (t0, t1, N, A, [1, 0])
    us, timesteps, h = integrate(prob)
    positions = us[1, :]
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

let
counter = 0
for γ in (0.5, 1, 2, 4) # Un, under, critical, over
    for timestep in (1, 1e-1, 1e-2, 1e-3)
        counter += 1
        timesteps, pos_est, pos_ana, resids, errors = task_2d(γ, timestep)
        fig, ax, plt = lines(timesteps, pos_est, label="Estimated position")
        ax.title = "γ = $γ, h = $timestep"
        lines!(timesteps, pos_ana, label="Analytical position")
        lines!(timesteps, resids, label="Residuals")
        lines!(timesteps, errors, label="Errors")
        axislegend()
        fig |> display
        @info "Finnished with $counter/16 plots"
    end
end
end