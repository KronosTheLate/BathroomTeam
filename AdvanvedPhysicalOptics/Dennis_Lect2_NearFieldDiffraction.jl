using GLMakie; Makie.inline!(true)

λ = 633e-9  # Wavelength
k = 2π/λ

x₁ᶥ_lim = -10e-3 # 4x₁ᶥ          # Lower calculation limit
x₂ᶥ_lim =  15e-3 # 4x₂ᶥ          # Upper calculation limit

D = 2                      # Distance from object to aperture plane
Dᶥ = 1.5                    # Distance from aperture plane to image plane
# d = 1/(1/D + 1/Dᶥ)        # Was defined in slides. Used in defining ηₓ
x̃₁ = -1e-3                  # Distance from optical axis to lower aperture edge
x̃₂ = -x̃₁                    # Distance from optical axis to upper aperture edge
x₁ᶥ = x̃₁/D * (D+Dᶥ)         # Projection of x̃₁ᶥ in image plane  x̃₁/D == x₁ᶥ/(D+Dᶥ)  as they are tan of same angle
x₂ᶥ = x̃₂/D * (D+Dᶥ)         # Projection of x̃₂ᶥ in image plane

F = √(/(λ*Dᶥ*(D+Dᶥ), 2D))
a₁(xᶥ) = (x₁ᶥ-xᶥ)/F         # Lower integration limit for ηₓ
a₂(xᶥ) = (x₂ᶥ-xᶥ)/F         # Upper integration limit for ηₓ
N = abs(x̃₁*x̃₂)/(λ*D)        # From first slide, N=x̃^2/(λ*D) - Fresnel zones?
@show N

# x̃ᵤ(xᶥ) = xᶥ/(D+Dᶥ)*D            # Intersect of line between source point (at optical axis) and aperture plane
# ηₓ(xᶥ) = √(2/(λ*d))*(x)

function integrate(f, lower, upper; step)
    argument_values = range(lower, upper; step)
    return sum(f, argument_values) * step
end

# using Integrals
# let 
#     prob = IntegralProblem((x, p) -> identity(x), 0, 100)
#     sol = solve(prob, QuadGKJL())
# end


f(ηₓ) = cis(-π/2 * ηₓ^2)
E(xᶥ) = √(im/2)*integrate(f, a₁(xᶥ), a₂(xᶥ), step= F/x₂ᶥ / 100)   # step=δη << F/x₂ᶥ
S(xᶥ) = abs2(E(xᶥ))


let 
    fig = Figure()
    ax = Axis(fig[1, 1])
    ax.xlabel = "z"
    ax.ylabel = "x"
    scatter!([0], [0], label="Object", markersize=30)
    begin # Plotting propegation circles around object. Did not work - axes are not of equal scale
        # ϕ = range(-π/2, π/2, 50)
        # total_y_height = x₂ᶥ_lim - x₁ᶥ_lim
        # circ_attenuator = 0.1total_y_height
        # lines!(cos.(ϕ).*circ_attenuator, sin.(ϕ) .* circ_attenuator, label="Test")
    end
    total_y_height = x₂ᶥ_lim - x₁ᶥ_lim
    lines!([D, D], [x₁ᶥ_lim, x̃₁], label="Aperture")
    lines!([D, D], [x̃₂, x₂ᶥ_lim], label="Aperture", color=Cycled(1))
    vlines!([D+Dᶥ], label="Image plane", linestyle="-", color=:black)


    N_points = 1000
    
    display_with_z = (D+Dᶥ)

    xᶥs = range(x₁ᶥ_lim, x₂ᶥ_lim, N_points)
    intensities = S.(xᶥs)
    intensities_normalized = intensities ./ maximum(intensities)
    intensity_coordinates = @. D + Dᶥ + display_with_z*intensities_normalized
    lines!(intensity_coordinates, xᶥs, label="Intensity")


    # vlines!([D], ymin=0, ymax=0.5)
    axislegend(position=(0, 1), merge=true)
    ylims!(x₁ᶥ_lim, x₂ᶥ_lim)
    fig |> display
end


##

let 
    
    N_points = 200
    
    xᶥs = range(x₁ᶥ_lim, x₂ᶥ_lim, N_points)
    intensities = S.(xᶥs)
    fig, ax, plt = scatterlines(xᶥs, intensities)
    fig |> display

end

#¤ For small enough slit, do we get sinc²?