#¤ Class problem 1: Linear eigenproblems

using LinearAlgebra
A = [
    1 2 3 4
    1 5 6 7
    1 1 8 9
    1 1 1 0
]

##

function fixed_point_iteration(A::AbstractMatrix; atol=1e-8, maxiters::Integer = 100)
    counter = 0

    #? Input check
    isequal(size(A)...)  ||  throw(ArgumentError("Input matrix A is not square.\nIt has size $(size(A))"))

    #? 1) Pick arbitrary vector:
    y = rand(size(A)[1])
    while true
        if counter ≥ maxiters
            @warn "Reached maximal number of iterations ($maxiters).\nTerminating."
            break
        end
        ##? 2) Calculate x = Ay and normalize x
        x = A*y
        x = x / norm(x)

        #? 3) If the change ||x-y|| is sufficiently small, terminate
        Δ = norm(x-y)
        println("Δ = $(round(Δ, sigdigits=10))")
        if Δ ≤ atol
            @info "Δ ≤ atol = $atol\nTerminating after $counter iterations."
            break
        end

        #? 4) Otherwise set y = x and go back to 2)
        y = x
        counter += 1
    end
    λ = (y⋅(A*y)) / (y⋅y)
    return (vec=x, val=λ)
end

function inverse_iteration(A::AbstractMatrix; atol=1e-8, maxiters::Integer = 100)
    local A⁻¹
    try
        A⁻¹ = inv(A)
    catch e
        @warn "Encountered error in inverting A.\nThrowing error."
        rethrow(e)
    end

    vec, val = fixed_point_iteration(A⁻¹; atol, maxiters)
    val = inv(val)
    return (;vec, val)
end

function shift_and_invert(A::AbstractMatrix, μ::Number; atol=1e-8, maxiters::Integer = 100)
    local Ã
    try
        Ã = inv(A-μ*LinearAlgebra.I)
    catch e
        @warn "Encountered error in inverting A-μ.\nThrowing error."
        rethrow(e)
    end

    vec, val = fixed_point_iteration(Ã; atol, maxiters)
    val = μ + inv(val)
    return (;vec, val)
end

inverse_iteration(A)
fixed_point_iteration(A)
shift_and_invert(A, -1.5)

##

B = 1/5 * [
    3 4 0
   -4 3 0
    0 0 5
]

fixed_point_iteration(B)
inverse_iteration(B)
shift_and_invert(B, -1e10)
shift_and_invert(B, 100000000000)
shift_and_invert(B, 10*cis(deg2rad(180)))


##¤ Class problem 2: Nonlinear problems
#? Defining the nonlinear function:

#* The problem: f(x⃗)x⃗ = b⃗

function relaxation_method(f, b⃗, x⃗₀; deflation=0, atol=1e-8, maxiters::Integer=100)
    0 ≤ deflation ≤ 1  ||  throw(ArgumentError("The deflation parameter is not between 0 and 1 (it is $deflation)."))
    counter = 0
    x⃗ = x⃗₀  # init
    while true
        if counter ≥ maxiters
            @warn "Reached maximal number of iterations ($maxiters).\nTerminating."
            break
        end
        x⃗_old = x⃗
        f_evaluated = f(x⃗_old)
        # x⃗ = f_evaluated \ b⃗
        y⃗ = f_evaluated \ b⃗
        x⃗ = (1-deflation)*y⃗ + deflation * x⃗_old
        counter += 1
        Δ = norm(x⃗ - x⃗_old)
        println("||x⃗ - x⃗_old|| = $(round(Δ, sigdigits=10))")
        if Δ ≤ atol
            @info "Δ ≤ atol = $atol\nTerminating after $counter iterations."
            break
        end
    end
    return x⃗
end

let 
    f(x⃗) = [1 2; 3 4] + norm(x⃗)*[0 1; -1 0]
    b⃗ = [0, 1/2]
    x⃗₀ = [0, 0]
    relaxation_method(f, b⃗, x⃗₀, deflation=0.7)
end

let 
    f(x⃗) = [1 2; 3 4] + norm(x⃗)*[0 1; -1 0]
    b⃗ = [0, 0.66]
    x⃗₀ = [0, 0]
    relaxation_method(f, b⃗, x⃗₀, deflation=0.91, maxiters=10^3)
end
