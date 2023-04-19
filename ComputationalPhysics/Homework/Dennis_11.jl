if occursin("dennishb", homedir())
    using Pkg
    Pkg.activate("ComputationalPhysics", shared=true)
end
using Symbolics
using GLMakie; Makie.inline!(true)

𝒢_Φ(x, xᶥ, y, yᶥ) = -1/4π * log(Complex((x-xᶥ)^2 + (y-yᶥ)^2))
let 
    @variables x xᶥ y yᶥ
    Dx = Differential(x)
    Dy = Differential(y)
    𝒢_Φ_expr = 𝒢_Φ(x, xᶥ, y, yᶥ)

    𝒢_Ex_expr = -Dx(𝒢_Φ_expr) |> expand_derivatives
    𝒢_Ex = build_function(𝒢_Ex_expr, [x xᶥ y yᶥ])|>eval
    𝒢_Ex(x, xᶥ, y, yᶥ) = 𝒢_Ex([x, xᶥ, y, yᶥ])

    𝒢_Ey_expr = -Dy(𝒢_Φ_expr) |> expand_derivatives
    𝒢_Ey = build_function(𝒢_Ey_expr, [x xᶥ y yᶥ])|>eval
    𝒢_Ey(x, xᶥ, y, yᶥ) = 𝒢_Ey([x, xᶥ, y, yᶥ])
    # blay = -Dy(greenfunc) |> expand_derivatives
    # blax, blay
end
#=
let 
    @variables x y
    testexpr = x^2 + 2y
    f_expr = build_function(testexpr, [x, y])|>eval
    f_expr(a, b) = f_expr([a, b])
    f_expr([3, 1]), f_expr(3, 1)
end
=#
##

let #¤ a1) Plot the potential in the domain
    N = 1000
    xmin, xmax = -2, 2
    xs = range(xmin, xmax, N)
    ymin, ymax = -2, 2
    ys = range(ymin, ymax, N)
    xᶥ = 1
    yᶥ = 0
    𝒢s = [𝒢_Φ(x, xᶥ, y, yᶥ) for x in xs, y in ys]
    heatmap(xs, ys, abs.(𝒢s))
end
function parameters(N)
    a=1
    b=1
    return (;N, a, b)
end

E_ext(x, y, p) = p.a*x + p.b*y
