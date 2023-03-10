##¤ Preamble
let  # I have started using enviroments. Do not mind this if you dont.
    using Pkg
    active_project = splitpath(Pkg.project().path)[end-1]
    expected_project = "ComputationalPhysics"
    active_project != expected_project  &&  @warn "Active project is not $expected_project"
end
using DifferentialEquations
using GLMakie; Makie.inline!(false)
using LinearAlgebra  # For norm, normalize


##¤ Solving 1a)
#* u = r
function ddu_particleincell!(u̇̇, u̇, u, p, t)
    m = p.m
    R = p.R(t)
    r = u
    u̇̇ .= 1/m * (-norm(r - R, 2)) * (r - R)
end
u0 = rand(3) .- 0.5
u̇0 = rand(3) .- 0.5
tspan = (0.0, 100.0)
p = (R=t->t≤25 ? [0, 0, 0] : [1, 1, 1], m=1)
prob_partitioned = SecondOrderODEProblem(ddu_particleincell!, u̇0, u0, tspan, p)
##
N_steps = 10^3
Δt = (tspan[end]-tspan[begin])/N_steps
sol = solve(prob_partitioned, VelocityVerlet(), dt=Δt)
rx = sol[1, :]
ry = sol[2, :]
rz = sol[3, :]
plt = lines(rx, ry, rz, label="r⃗", axis=(xlabel="x", ylabel="y", zlabel="z"))

##
fig = Figure()
sg = SliderGrid(fig[2, 1], (label="ind", range=eachindex(sol.t), ))
ind = sg.sliders[1].value
ax = Axis(fig[1, 1], title=@lift("t = $(round(sol.t[$ind], digits=2))"))
limits!(ax, extrema(rx), extrema(ry))
point = @lift Point2f.(rx[$ind[]], ry[$ind[]])#, rz[$ind[]])
plt = scatter!(point)
fig

##
# ind_step = (tspan[end]-tspan[begin])/n_frames / Δt |> round |> Int
duration = 5
fps = 100
n_frames = fps * duration
for i in range(1, length(sol.t), n_frames)
    ind[] = round(Int64, i)
    sleep(1/fps)
end



##¤ Implementing Kepler problem with DifferentialEquations.jl
#* u:   State vector.
#* p:   Parameters of simulation.
#* t:   Timepoint.
function du_kepler(u, p, t)
    length(u) == 4  ||  @warn "State vector is not of length 4, as expected."
    rx, ry, px, py = u
    pos_mag³ = hypot(rx, ry)^3
    return [px, py, -rx/pos_mag³, -ry/pos_mag³]
end

function ddu_kepler(uᶥ, u, p,t)  # In this form, u is r (position vector)
    rx, ry = u
    pos_mag³ = hypot(rx, ry)^3
    uᶥᶥ = [-rx/pos_mag³, -ry/pos_mag³]
    return uᶥᶥ
end

r0 = [1.0, 0.0]
p0 = [0, 1/√r0[1]]
u0 = vcat(r0, p0)
tspan = (0.0, 100.0)
prob = ODEProblem(du_kepler, u0, tspan)
uᶥ0 = p0
prob_partitioned = SecondOrderODEProblem(ddu_kepler, uᶥ0, u0[1:2], tspan)

##


N_steps = 10^3
Δt = (tspan[end]-tspan[begin])/N_steps

sol_euler = solve(prob, Euler(), dt=Δt)
sol_velocityverlet = solve(prob_partitioned, VelocityVerlet(), dt=Δt)

u_ref(t, Ω=1) = [cos(Ω*t), sin(Ω*t), -sin(Ω*t), cos(Ω*t)]
sol_ref = hcat(u_ref.(sol_euler.t)...)
rx_ref = sol_ref[1, :]
ry_ref = sol_ref[2, :]
px_ref = sol_ref[3, :]
py_ref = sol_ref[4, :]
r_mag_ref = hypot.(rx_ref, ry_ref)

residuals = sol_euler .- sol_ref
rx_resid = residuals[1, :]
ry_resid = residuals[2, :]

# sol = sol_euler
sol = sol_velocityverlet

alg = sol.alg
alg_string = string(alg)[begin:end-2]
plot_title = "Algorithm = $alg_string\n$N_steps steps"
rx = sol[1, :]
ry = sol[2, :]
px = sol[3, :]
py = sol[4, :]
r_mag = hypot.(rx, ry)
markersizes = sol.t / maximum(sol.t)  .* 2 .+ 0.5

let 
    plt1 = plot(sol.t, r_mag, label="|r⃗|", xlabel="t")
    plot!(sol.t, r_mag_ref, label="|r⃗|_ref")
    plt2 = plot(rx, ry, label="r⃗", lw=markersizes, xlabel="r_x", ylabel="r_y")
    plot!(rx_ref, ry_ref, label="r⃗_ref", lw=markersizes, aspect_ratio=:equal)
    plot(plt1, plt2, layout=(1, 2); plot_title, plot_titlevspan=0.15)|>display
end

##
for N_steps in [100*2^i for i in 0:10]
    Δt = (tspan[end]-tspan[begin])/N_steps

    sol_euler = solve(prob, Euler(), dt=Δt)
    sol_velocityverlet = solve(prob_partitioned, VelocityVerlet(), dt=Δt)

end

##¤ ModelingToolkit:
using DifferentialEquations
using ModelingToolkit
@variables t rx(t) ry(t)
∂ₜ = Differential(t)
∂ₜ² = ∂ₜ^2

@named kepler_eqns = ODESystem([∂ₜ²(rx) .~ -rx/hypot(rx, ry)^3,
                                ∂ₜ²(ry) .~ -ry/hypot(rx, ry)^3]
)
kepler_eqns_simplified = structural_simplify(kepler_eqns)
u0 = [rx=>1, ry=>0, ∂ₜ(rx)=>0, ∂ₜ(ry)=>1]
tspan = (0.0, 10.0)
prob = ODEProblem(kepler_eqns_simplified, u0, tspan)
prob = SecondOrderODEProblem(kepler_eqns_simplified, u0, tspan)
prob = DynamicalODEProblem(kepler_eqns, u0, tspan)
prob = HamiltonianProblem(kepler_eqns, u0, tspan)
sol = solve(prob, Tsit5())
sol.u