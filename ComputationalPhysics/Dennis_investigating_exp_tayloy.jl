using TaylorSeries
using GLMakie; Makie.inline!(true)

let
    fig = Figure()
    ax = Axis(fig[1, 1])
    ts = range(0, 40, 100)
    for n in 5:10:90
        t = Taylor1(n)
        f = exp(t)
        lines!(ts, f.(ts), label="$n terms")
    end
    lines!(ts, exp, label="exp(x)", color=:black)
    axislegend(position=(0, 1))
    ax.yscale=log
    fig |> display
end