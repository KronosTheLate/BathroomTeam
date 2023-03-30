using Symbolics
@variables x B
Dx = Differential(x)

once_differentiated = Dx(1/cosh(B*x)) |> expand_derivatives |> simplify
Dx(once_differentiated) |> expand_derivatives |> simplify

##¤

using DifferentialEquations
using FFTW
using GLMakie; Makie.inline!(true); update_theme!(resolution=(1000, 1000), fontsize=30)
using LinearAlgebra     # For norm

function get_lims(data...)  # Useful for setting common colorbars
    extremas = extrema.(data)
    lowest = minimum(getindex.(extremas, 1))
    highest = maximum(getindex.(extremas, 2))
    return (lowest, highest)
end

"""
k₀ = 2π/x_span
"""
function approximate_derivative(evaluated_points::AbstractVector, k₀::Real; tol=1e-10, check=false)
    dft = fft(evaluated_points)
    N = length(evaluated_points)
    prefactors = fftfreq(N, N) .* im*k₀
    output = ifft(prefactors .* dft)
    if check
        max_imag_part = abs(maximum(imag.(output)))
        max_imag_part ≥ tol  &&  @warn("Largest imaginary part of derivative is $(round(max_imag_part, sigdigits=5))")
    end
    return output
end

"""
k₀ = 2π/x_span
"""
function approximate_2nd_derivative(evaluated_points::AbstractVector, k₀::Real; tol=1e-10, check=false)
    dft = fft(evaluated_points)
    N = length(evaluated_points)
    prefactors = (fftfreq(N, N) .* im*k₀).^2
    output = ifft(prefactors .* dft)
    if check
        max_imag_part = abs(maximum(imag.(output)))
        max_imag_part ≥ tol  &&  @warn("Largest imaginary part of derivative is $(round(max_imag_part, sigdigits=5))")
    end
    return output
end


##¤ Function to get a NamedTuple of parameters. `check` and `tol` are for approximate_derivative
function parameters(N::Integer; A=1, tol=2, check=false)  # Function, because want to update for convergence study
    L = 50
    x_min = -L/2
    x_max = L/2
    xs = range(x_min, x_max, N+1); xs=xs[begin:end-1]
    t_max = 20π
    tspan = (0, t_max)
    dt = 1/N * 1.0
    k₀ = 2π/(x_max-x_min)

    B = A/√2
    Ω = B^2
    return (;tol, check, A, B, Ω, xs, k₀, tspan, dt)
end

u_ref_func(t, p) = @. p.A/cosh(p.B * p.xs) * cis(-p.Ω*t)
post_processor(solution) = reduce(hcat, solution) .|> abs2 |> transpose
log10_offset(x) = log10(x+1e-7)

function update_slope_schr_full!(du, u, p, t)  #¤ For task b
    ∂ₓ²u = approximate_2nd_derivative(u, p.k₀; tol=p.tol, check=p.check)
    @. du = 1/im * (∂ₓ²u + abs2(u)*u)
    return nothing
end
function update_slope_schr_dispersive!(du, u, p, t)  #¤ For task c
    ∂ₓ²u = approximate_2nd_derivative(u, p.k₀; tol=p.tol, check=p.check)
    @. du = 1/im * (∂ₓ²u)
    return nothing
end
function update_slope_schr_NL!(du, u, p, t)  #¤ For task c
    @. du = 1/im * (abs2(u)*u)
    return nothing
end

##¤ b) Building on last week's script, solve new problem for L=50, N=200, tₘₐₓ = 20π
p = parameters(200, A=√2)
u0 = u_ref_func(0, p)
prob = ODEProblem(update_slope_schr_full!, u0, p.tspan, p)
sol = solve(prob, RK4(); adaptive=false, p.dt)
u_ref = [u_ref_func(t, p) for t in sol.t]
u_resid = sol.u - u_ref 

# scale = identity
scale = log10
# scale = log10_offset

sol_u_processed = sol.u |> post_processor .|> scale
u_ref_processed = u_ref |> post_processor .|> scale
u_resid_processed = u_resid |> post_processor .|> scale

let  #¤ Plotting
    clims = get_lims(sol_u_processed, u_ref_processed)
    fig, ax, hm1 = heatmap(sol.t, p.xs, sol_u_processed, colorrange=clims)
    ax.ylabel="x"
    ax.xlabel="t"
    ax.title = "Numerical solution"
    ax2, hm2 = heatmap(fig[2, 1], sol.t, p.xs, u_ref_processed, colorrange=clims)
    ax2.ylabel="x"
    ax2.xlabel="t"
    ax2.title = "Reference solution"
    ax3, hm3 = heatmap(fig[3, 1], sol.t, p.xs, u_resid_processed)#, colorrange=clims)
    ax3.ylabel="x"
    ax3.xlabel="t"
    ax3.title = "Residuals"
    Colorbar(fig[1:2, 2], colorrange=clims)
    Colorbar(fig[3, 2], hm3)
    Label(fig[0, :], string("Scale = ", scale))
    display(fig)
end

##¤ c) Remove nonlinear term |a|²a, use gaussian initial conditions

p = parameters(200, A=√2)
u0 = u_ref_func(0, p)

#¤ Gaussians
# u0 = ComplexF64[exp(-0.25x^2) for x in p.xs]
# u0 = ComplexF64[exp(-0.25(x-4)^2) for x in p.xs]
# u0 = ComplexF64[exp(-0.5(x-4)^2) for x in p.xs]   #¤ Not even changing shape, just oscillatin
# u0 = ComplexF64[exp(-2x^2) for x in p.xs]         #¤ "Breathing"/"pumping"
# u0 = ComplexF64[exp(  -2(x-4)^2) for x in p.xs]   #¤ "Pumping" + oscillation
# u0 = ComplexF64[exp(  -2(x-4)^2) - exp(  -2(x+4)^2) for x in p.xs]
# u0 = ComplexF64[exp(  -2(x-4)^2) + exp(  -2(x+4)^2) for x in p.xs]

#¤ From back of sheet
# u0 = ComplexF64[10*exp(-1/30 * (x-5)^2) for x in p.xs]
# u0 = @. √8 / cosh(2*p.xs)                       |> ComplexF64
# u0 = @. √2*exp(0.1im*p.xs) / cosh(p.xs)         |> ComplexF64
u0 = @. √8*exp(-0.2im*p.xs) / cosh(2*p.xs+20) + √8*exp(0.1im*p.xs) / cosh(2*p.xs-30) - √8/cosh(2*p.xs)   |> ComplexF64
# u0 = @. √8 / cosh(p.xs)                         |> ComplexF64


scale = log10_offset
for update_slope! in (update_slope_schr_full!, update_slope_schr_dispersive!, update_slope_schr_NL!)|>reverse
    prob = ODEProblem(update_slope!, u0, p.tspan, p)
    sol = solve(prob, RK4(); adaptive=false, p.dt)
    u_ref = [u_ref_func(t, p) for t in sol.t]
    u_resid = sol.u - u_ref 
    sol_u_processed = sol.u |> post_processor .|> scale
    u_ref_processed = u_ref |> post_processor .|> scale
    u_resid_processed = u_resid |> post_processor .|> scale

    fig, ax, hm = heatmap(sol.t, p.xs, sol_u_processed)
    ax.ylabel="x"
    ax.xlabel="t"
    ax.title = "Magnitude"
    ax2, hm2 = heatmap(fig[2, 1], sol.t, p.xs, reduce(hcat, sol.u) |>transpose .|> angle .|> rad2deg)
    ax2.ylabel="x"
    ax2.xlabel="t"
    ax2.title = "Phase"
    Label(fig[0, :], string("Scale = ", scale, "\nEquation = ", update_slope!), tellwidth=false)
    Colorbar(fig[1, 2], hm)
    Colorbar(fig[2, 2], hm2)
    display(fig)
end












##¤=== Stuff fron last week ===¤##
##¤ Convergence study

function L∞∞(N)
    p = parameters(N; tol=100)
    (;xs, tspan, dt) = p
    u0 = u_ref_func(0, xs)
    prob_convergence = ODEProblem(update_slope_schr!, u0, tspan, p)
    sol = solve(prob_convergence, RK4(); adaptive=false, dt)
    u_ref = [u_ref_func(t, xs) for t in sol.t]
    u_resid = sol.u - u_ref 

    return norm(norm.(u_resid, Inf), Inf)
end

function L∫∫(N)
    p = parameters(N; tol=100)
    (;xs, tspan, dt) = p
    u0 = u_ref_func(0, xs)
    prob_convergence = ODEProblem(update_slope_schr!, u0, tspan, p)
    sol = solve(prob_convergence, RK4(); adaptive=false, dt)
    u_ref = [u_ref_func(t, xs) for t in sol.t]
    u_resid = sol.u - u_ref 

    h = step(p.xs)  # Integrating all time and space
    integrated_over_space = [h * sum(abs.(residuals)) for residuals in u_resid]
    integrated_over_time = dt * sum(abs.(residuals) for residuals in integrated_over_space)
    return abs(integrated_over_time)
end
##
# Ns = [round(Int, 100*1.2^i) for i in 0:5]
Ns = 20:2:150
dts = getproperty.(parameters.(Ns), :dt)
errorfunc = L∞∞
# errorfunc = L∫∫

# scale = identity
# scale = log10
# scale = log10_offset

errors = errorfunc.(Ns)
mask_convergence = isfinite.(errors)
let
    fig, ax, plt = lines(Ns[mask_convergence], errors[mask_convergence])
    ax.xlabel = "N"
    ax.ylabel = "Error measure"
    ax.yscale = log10
    ax.xscale = log10
    Label(fig[0, :], string("Convergence study\n","Error function = ", errorfunc), tellwidth=false)
    display(fig)
end
# energy_expectationval = [ψ' * Hψ(ψ, p) for ψ in sol.u]


##¤ d) Different initial conditions

p = parameters(100)
length(p.xs)
u0 = u_ref_func(0, p.xs)
# u0 = ComplexF64[exp(-0.25x^2) for x in p.xs]
# u0 = ComplexF64[exp(-0.25(x-4)^2) for x in p.xs]
# u0 = ComplexF64[exp(-0.5(x-4)^2) for x in p.xs]   #¤ Not even changing shape, just oscillatin
# u0 = ComplexF64[exp(-2x^2) for x in p.xs]         #¤ "Breathing"/"pumping"
# u0 = ComplexF64[exp(  -2(x-4)^2) for x in p.xs]   #¤ "Pumping" + oscillation
# u0 = ComplexF64[exp(  -2(x-4)^2) - exp(  -2(x+4)^2) for x in p.xs]
# u0 = ComplexF64[exp(  -2(x-4)^2) + exp(  -2(x+4)^2) for x in p.xs]
prob = ODEProblem(update_slope_schr!, u0, p.tspan, p)
sol = solve(prob, RK4(); adaptive=false, p.dt)

# scale = log10
scale = log10_offset
sol_u_processed = sol.u |> post_processor .|> scale
let  #¤ Plotting
    fig, ax, hm = heatmap(p.xs, sol.t, sol_u_processed)
    ax.xlabel="x"
    ax.ylabel="t"
    ax.title = "Numerical solution"
    ax2, hm2 = heatmap(fig[2, 1], p.xs, sol.t, angle.(reduce(hcat, sol.u)))
    Label(fig[0, :], string("Scale = ", scale), tellwidth=false)
    Colorbar(fig[1, 2], hm)
    Colorbar(fig[2, 2], hm2)
    display(fig)
end