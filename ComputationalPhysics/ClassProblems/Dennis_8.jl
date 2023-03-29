#?  Classroom problem: Differential operator in Fourier basis
#?  Consider the interval [0, L] with periodic boundary conditions and
#?  construct a MATLAB function handle that performs the derivative on
#?  any function sampled in this interval:

#¤  1. Discretize the interval into N points, please note that the points
#¤  x = 0 and x = L are the same point in a periodic system and should
#¤  not appear twice.

#*  2. Approximate the derivative by first switching to Fourier space. Then
#*  multiply each Fourier component with the respective −ik. Finally,
#*  switch back to real space.

#?  Keep the length L as a parameter.

#?  Select L = 2π and test your implementation for different values of N and
#?  a number of functions, for example
#?  f(x) = cos(2x)
#?  f(x) = cos(1.7x)
#?  f(x) = x
#?  f(x) = x²
#?  f(x) = exp[−(x − π)² ]
#?  f(x) = exp[−(x − π)² /4]
#?  Also find the analytic derivatives and compare.

##¤ 1
N = 100     #* Number of points
L = 2π
xs = range(0, length=N, step=L/N)    # End exlusive range from 0 to L with N points
k₀ = 2π/L
##* 2
using FFTW
using ForwardDiff
using GLMakie; Makie.inline!(true)

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

# approximate_derivative([sin(x) for x in range(0, 2π, 101)[begin:end-1]], 1)|>display
approximate_derivative(cospi.(2*xs), k₀)|>display

function plot_true_and_approx_derivs(f, xs=xs, k₀=k₀)
    evaluated_points = f.(xs)
    approx_deriv = approximate_derivative(evaluated_points, k₀) .|> real
    true_deriv = ForwardDiff.derivative.(f, xs)
    err = approx_deriv .- true_deriv
    fig, ax, _ = lines(xs, f.(xs), label="f(x)")
    ax.xlabel = "x"
    lines!(xs, approx_deriv, label="Approx ∂ₓf")
    lines!(xs, true_deriv, label="True ∂ₓf", linestyle=:dash)
    approx_2nd_deriv = approximate_2nd_derivative(evaluated_points, k₀) .|> real
    lines!(xs, approx_2nd_deriv, label="Approx ∂ₓ²")
    true_2nd_deriv = ForwardDiff.derivative.(vals -> ForwardDiff.derivative(f, vals), xs)
    lines!(xs, true_2nd_deriv, label="True ∂ₓ²", linestyle=:dash)

    axislegend()
    display(fig)
    return fig, ax
end
plot_true_and_approx_derivs(x->cos(2x))
plot_true_and_approx_derivs(x->sin(x))
plot_true_and_approx_derivs(identity)
ylims!(current_axis(), -5, 5)
display(current_figure())
