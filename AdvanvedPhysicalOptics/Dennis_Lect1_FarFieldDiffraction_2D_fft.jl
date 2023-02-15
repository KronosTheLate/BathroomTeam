##¤ Lecture 2 - write program that computes diffraction pattern of rectangular slit
using GLMakie; Makie.inline!(true); update_theme!(colormap=:greys)

λ = 550e-9
k = 2π/λ

begin   #? Defining slit by transmission function τ̃.
        #? Note that τ̃ gives fraction of maximal power transmitted.
    slit_Δx::Float64 = 10λ        # Slit is centered
    slit_Δy::Float64 = 1λ      # on x=0, y=0
    τ̃_rectangular(x̃, ỹ) = ifelse((abs(x̃) ≤ slit_Δx/2) && (abs(ỹ) ≤ slit_Δy/2), 1 , 0)
    τ̃ = τ̃_rectangular
end

N::Int64 = 100   # Number of points in x and y dimension, both angle and space.
input_sidelength = 1.5*max(slit_Δx, slit_Δy)
x̃s = range(-input_sidelength/2, input_sidelength/2, N)
ỹs = range(-input_sidelength/2, input_sidelength/2, N)

begin #? Visualizing the transmission function τ̃
    transmission_matrix = [τ̃(x̃, ỹ) for x̃ in x̃s, ỹ in ỹs]
    heatmap(x̃s, ỹs, transmission_matrix, 
        axis=(xlabel="x", ylabel="y", title="Slit",
         xticks = MultiplesTicks(6, λ, "λ"), yticks = MultiplesTicks(6, λ, "λ"))
    ) |> display
end

function my_2D_fft(α, β; x̃s=x̃s, ỹs=ỹs, A::Float64=1.0, k::Float64=k, τ̃=τ̃_rectangular)::ComplexF64 # A = Amplitude, E_in/R in final slide of Lecture 2
    A_normalized = A / length(x̃s) / length(ỹs)
    output_nonnormalized = sum(τ̃(x̃, ỹ) * cis(k*(α*x̃ + β*ỹ)) for x̃ in x̃s, ỹ in ỹs)
    return A_normalized * output_nonnormalized
end

αs = range(-deg2rad(50), deg2rad(50), N)
βs = range(-deg2rad(50), deg2rad(50), N)
@time response = [my_2D_fft(α, β) for α in αs, β in βs]

let #? Visualizing the diffraction pattern
    fig, ax1, plt1 = heatmap(x̃s, ỹs, zs, 
        axis=(xlabel="x", ylabel="y", title="Slit",
         xticks = MultiplesTicks(6, λ, "λ"), yticks = MultiplesTicks(6, λ, "λ"))
    )
    amplitudes = abs.(response)
    heatmap(fig[1, 2], αs, βs, amplitudes,
        axis=(xlabel="α", ylabel="β", title="Diffraction pattern",
        #  xticks = MultiplesTicks(6, λ, "λ"), yticks = MultiplesTicks(6, λ, "λ")
    ))
    fig |> display
end

# u = -α/λ
# v = -β/λ

##? using FFTW
using FFTW

let #? Visualizing the diffraction pattern
    fig, ax1, plt1 = heatmap(x̃s, ỹs, zs, 
        axis=(xlabel="x", ylabel="y", title="Slit",
         xticks = MultiplesTicks(6, λ, "λ"), yticks = MultiplesTicks(6, λ, "λ"))
    )
    amplitudes = abs.(fftshift(fft(transmission_matrix)))

    heatmap(fig[1, 2], fftshift(fftfreq(length(x̃s)), 1/step(x̃s)), fftshift(fftfreq(length(ỹs), 1/step(ỹs))), amplitudes,
        axis=(xlabel="α", ylabel="β", title="Diffraction pattern",
        #  xticks = MultiplesTicks(6, λ, "λ"), yticks = MultiplesTicks(6, λ, "λ")
    ))
    fig |> display
end
