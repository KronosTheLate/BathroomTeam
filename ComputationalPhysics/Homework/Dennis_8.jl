using DifferentialEquations
using FFTW  #* Function from class:
using GLMakie; Makie.inline!(true)
using ForwardDiff

##
"""
k₀ = 2π/x_span
"""
function approximate_derivative(evaluated_points::AbstractVector, k₀::Real; tol=1e-10)
    dft = fft(evaluated_points)
    N = length(evaluated_points)
    prefactors = fftfreq(N, N) .* im*k₀
    output = ifft(prefactors .* dft)
    max_imag_part = abs(maximum(imag.(output)))
    max_imag_part ≥ tol  &&  @warn("Largest imaginary part of derivative is $(round(max_imag_part, sigdigits=5))")
    return output
end

function approximate_2nd_derivative(evaluated_points::AbstractVector, k₀::Real; tol=1e-10)
    dft = fft(evaluated_points)
    N = length(evaluated_points)
    prefactors = (fftfreq(N, N) .* im*k₀).^2
    output = ifft(prefactors .* dft)
    max_imag_part = abs(maximum(imag.(output)))
    max_imag_part ≥ tol  &&  @warn("Largest imaginary part of derivative is $(round(max_imag_part, sigdigits=5))")
    return output
end


##¤ Setup grid and input for `solve`
x_max = 10
x_min = -x_max
N = 101
xs = range(x_min, x_max, N)
t_max = 4π
tspan = (0, t_max)
dt = 1/(N-1)
k₀ = 2π/(x_max-x_min)


#¤ du: slope vector (will be updated)
#¤  u: current state
#¤  p: parameters
#¤  t: time dependence (not used)
function update_slope_schr!(du, u, p, t)
    (;tol, κ, ω, xs, k₀) = p
    k₀ = p.k₀
    xs = p.xs
    ∂ₓ²u = approximate_2nd_derivative(u, k₀; tol)
    @. du = 1/2im * (-∂ₓ²u + xs^2*u)
    return nothing
end
u_ref_func(t, xs) = @. exp(-p.κ*xs^2 - im*p.ω*t)
u0 = u_ref_func(0, xs)
p = (;tol=2, κ=0.5, ω=0.5, xs, k₀)
prob = ODEProblem(update_slope_schr!, u0, tspan, p)
sol = solve(prob, RK4())

post_processor(solution) = reduce(hcat, solution) .|> abs2# .|> log10
sol_u_processed = sol.u |> post_processor
u_ref = [u_ref_func(t, xs) for t in sol.t] |> post_processor
u_err = abs.(sol_u_processed .- u_ref)
function get_lims(data...)
    extremas = extrema.(data)
    lowest = minimum(getindex.(extremas, 1))
    highest = maximum(getindex.(extremas, 2))
    return (lowest, highest)
end

let  #¤ Plotting
    clims = get_lims(sol_u_processed, u_ref)
    fig, ax, hm1 = heatmap(xs, sol.t, sol_u_processed, colorrange=clims)
    ax.xlabel="x"
    ax.ylabel="time"
    ax.title = "Numerical solution"
    ax, hm2 = heatmap(fig[2, 1], xs, sol.t, u_ref, colorrange=clims)
    ax.xlabel="x"
    ax.ylabel="time"
    ax.title = "Reference solution"
    ax, hm3 = heatmap(fig[3, 1], xs, sol.t, u_err)#, colorrange=clims)
    ax.xlabel="x"
    ax.ylabel="time"
    ax.title = "Error"
    Colorbar(fig[1:2, 2], colorrange=clims)
    # Colorbar(fig[2, 2], hm2)
    Colorbar(fig[3, 2], hm3)
    #Colorbar(fig[:, 2], colorrange=clims)
    display(fig)
end