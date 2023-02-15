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
N = 50
γ = 0.01
A = [0 1; -1 -γ]
prob1 = (t0, t1, N, A, ones(size(A, 1)))
let
    us, timesteps, h = integrate(prob1)
    step_ = step(timesteps)
    us|>display
    # positions = first.(us)
    # velocities = last.(us)
    # fig, ax, plt = plot(timesteps, positions, label="Positions")
    # axislegend()
    # fig |> display
end

