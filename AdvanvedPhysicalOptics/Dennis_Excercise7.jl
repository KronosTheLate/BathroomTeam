if occursin("dennishb", homedir())
    using Pkg
    Pkg.activate("AdvancedPhysicalOptics", shared=true)
end
using GLMakie;
Makie.inline!(true);
update_theme!(resolution=(1920, 1080) .÷ 1.4, fontsize=40, linewidth=5, markersize=30)

##
using Unitful
let
    v = 4.2u"m/ms"
    f = 80u"MHz"
    Λ = v / f |> u"μm"

    λ = 632.8u"nm"
    L_max = 33u"mm"
    n = 2.26
    Q = /(2π * λ * L_max, n * Λ^2) |> NoUnits
end

##TODO  Diffraction efficiency for different sound wave amplitudes
using Interpolations
let
    input = [
        0.038
        0.091
        0.194
        0.29
        0.338
        0.404
        0.439
        0.486
        0.537
        0.589
        0.635
        0.69
        0.734
        0.788
        0.842
        0.891
        0.938
        0.99
        1.116
        1.246
    ]
    output = [
        0
        0
        0
        1.64
        3.7
        7.8
        10.3
        13.7
        17.5
        21.1
        24
        26.7
        28.6
        30.5
        31.7
        32.7
        33.3
        33.7
        34.6
        35.3
    ]
    global transfer_function = linear_interpolation(input, output, extrapolation_bc=Line())#, extrapolation_bc=Line())
    fig, ax, plt = scatter(input, output, label="Datapoints")
    ax.title = "Transfer function of AOM RF driver"
    ax.xlabel = "Driver Input  / V"
    ax.ylabel = "Driver Output  / V"
    xs = [-0.05; input; 1.5]
    lines!(xs, transfer_function.(xs), label="Linear interpolation + line extrapolation")
    axislegend(position=:rb)
    save(joinpath(homedir(), "Pictures", "transfer_curve.jpg"), fig)
    display(current_figure())
end

let  #* all units in Volts
    data = (
        driver_input=[
            0
            0.084
            0.198
            0.298
            0.366
            0.436
            0.485
            0.566
            0.635
            0.733
            0.83
            0.947
            1.097
            1.341
        ],
        order1=[
            0
            0
            0
            0
            0.023
            0.1
            0.175
            0.301
            0.45
            0.55
            0.63
            0.68
            0.74
            0.795
        ],
        order0=[
            1.025
            1.02
            0.97
            0.96
            0.94
            0.87
            0.85
            0.66
            0.55
            0.45
            0.36
            0.295
            0.25
            0.22
        ]
    )
    data = merge(data, (driver_output=transfer_function.(data.driver_input),))
    Δinput = 0.01
    Δoutput = 10e-3
    diffraction_efficiency = data.order1 ./ data.order0
    fig = Figure()
    ax = Axis(fig[1, 1], ylabel="Camera output / V", xlabel="Driver output  / V", rightspinevisible=false)
    scatterlines!(ax, data.driver_output, data.order1, label="∝ P₁", color=:black, marker=:rect)
    scatterlines!(ax, data.driver_output, data.order0, label="∝ P₀", color=:black, linestyle=:dash, marker=:circle)
    axislegend(position=:lc)

    ax2 = Axis(fig[1, 1], ylabel=rich("Diffraction efficiency ϵ", color=:red),
        yaxisposition=:right,
        yticklabelcolor=:red,
        rightspinecolor=:red,
        ytickcolor=:red)
    hidespines!(ax2, :l, :b, :t)

    scatterlines!(ax2, data.driver_output, diffraction_efficiency, label="Order 1 intensity / order 0 intensity", color=:red, marker=:utriangle)
    save(joinpath(homedir(), "Pictures", "diffraction_efficiency.jpg"), fig)
    display(fig)
end



let
    fig = Figure()
    ax = Axis(fig[1, 1])
    plt = lines!(ax, 0:10, 0:10)
    ylims!(ax, 1, 10)
    ax.yscale = log
    display(fig)
end

##TODO  Plot 0th and 1st order powers and the diffraction efficiency as a function of RF signal amplitude. [Dennis]

