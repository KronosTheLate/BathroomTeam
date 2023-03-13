dir_exc = "/home/dennishb/Uni/Semester/8. Sem/Advanced Physical Optics/Lab excercise 5"
using DelimitedFiles
using GLMakie; Makie.inline!(false)


r_s(n₁, n₂, θᵢ, θₜ) = /(n₁*cos(θᵢ) - n₂*cos(θₜ), 
                        n₁*cos(θᵢ) + n₂*cos(θₜ)
)
r_p(n₁, n₂, θᵢ, θₜ) = /(n₂*cos(θᵢ) - n₁*cos(θₜ), 
                        n₂*cos(θᵢ) + n₁*cos(θₜ)
)

δ(λ, n, d, θ) = 2π/λ * n*d*cos(θ)

##

let
    λ = 633e-9
    d1 = 100e-9

    n0 = 1
    n1 = 1.1
    n2 = 2
    θi_01 = 20|>deg2rad
    #¤  Snell's law:    n₁*sin(θᵢ) == n₂*sin(θₜ)
    θt_01 = n0*n1*sin(θi_01) |> asin; @info "01 went fine"
    θi_12 = θt_01
    θt_12 = n1*n2*sin(θi_12) |> asin; @info "01 went fine"
    r_01_p = r_p(n0, n1, θi_01, θt_01)
    r_12_p = r_p(n1, n2, θi_12, θt_12)
    r_01_s = r_s(n0, n1, θi_01, θt_01)
    r_12_s = r_s(n1, n2, θi_12, θt_12)

    δ1 = δ(λ, n1, d1, θt_01)

    rp = /(r_01_p + r_12_p*cis(2δ1), 1+r_01_p*r_12_p*cis(2δ1))
    rs = /(r_01_s + r_12_s*cis(2δ1), 1+r_01_s*r_12_s*cis(2δ1))
    @show rs
    @show rp
    rs_rp = rs/rp

    ψ(rs, rp) = atan(rs/rp)
    # rs_rp_ratio(ψ, Δ) = tan(ψ)*cis(Δ)

end


##
using Optim
xs = 1.0:10|>collect
ys = @. 5*xs^2 + rand()   # Coefficient 5

function lossfunc(parameters)
    α, β = parameters  # destructuring assignment. Makes next line more readable
    predictions = @. α*xs::Vector{Float64} + β*xs::Vector{Float64}^2
    return sum(abs, predictions-ys)
end
sol = optimize(lossfunc, [0, 0.0])
sol.minimizer
##
n_cauchy(A, B, C, λ) = A + B/λ^2 + C/λ^4
k_cauchy(D, E, F, λ) = D + E/λ^2 + F/λ^4

# form required by LsqFit
# λs is vector of indep vars from data
# p is vector of parameters that will be fitted
model_n(λs, p) = [n_cauchy(p[1], p[2], p[3], λ) for λ in λs]
ps = [633e-9]
dispersiondata_fused_silica = readdlm(joinpath(dir_exc, "dispersion_fused_silica.txt"))
dispersiondata_PMMA = readdlm(joinpath(dir_exc, "dispersion_PMMA.TXT"))
bla1, bla2 = eachcol(hcat(1:10, 11:20))
let
    data_p3 = readdlm(joinpath(dir_exc, "METALS_Gold_Johnson_Christy.txt"))
    λs, ns, ks = data_p3|>eachcol
    λs .*= 1e-6
    c = 3e8
    fs = c ./ λs
    ωs = 2π*fs

    ñ = ns .+ im .* ks
    ϵ_rs = @. ns^2 - ks^2
    ϵ_is = @. 2ns*ks
    ϵs = @. ϵ_rs + im*ϵ_is

    function lossfunc(p)
        ϵ∞, ωₚ, γ = p
        predictions_ϵ = [ϵ∞ - ωₚ^2/(ω^2+γ^2) + im*ωₚ^2*γ/(ω*(ω^2+γ^2)) for ω in ωs]
        residuals = predictions_ϵ - ϵs
        sum(abs2, residuals)
    end
    # sol = optimize(lossfunc, [1, 1e15, 1e14],)
    sol = optimize(lossfunc, [1, 1e15, 1e14], ParticleSwarm(n_particles=10))
    @show sol.minimum
    println("Minimizer")
    display(sol.minimizer)
    propertynames(sol)
end

(:method, :initial_x, :minimizer, :minimum, :iterations, :iteration_converged, :x_converged, :x_abstol, :x_reltol, :x_abschange, :x_relchange, :f_converged, :f_abstol, :f_reltol, :f_abschange, :f_relchange, :g_converged, :g_abstol, :g_residual, :f_increased, :trace, :f_calls, :g_calls, :h_calls, :ls_success, :time_limit, :time_run, :stopped_by)