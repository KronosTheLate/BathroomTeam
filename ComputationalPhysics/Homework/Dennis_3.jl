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
A
prob1 = (t0, t1, N, A, [1, 0])
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
    @show A
    @show B
    
    us_analytical = [A*exp(t*(-γ/2 + √(γ^2/4-1+0im))) + 
                        B*exp(t*(-γ/2 - √(γ^2/4-1+0im))) for t in timesteps
    ]
    lines!(timesteps, us_analytical.|>real, label="Analytical position")
    lines!(timesteps, velocities, label="Velocities")
    axislegend()
    # ax.xscale=log10
    fig |> display
end

[(√(x^2-1 + 0im) - x) / (2√(x^2-1 + 0im)) for x in range(00, 1000, 10000)]