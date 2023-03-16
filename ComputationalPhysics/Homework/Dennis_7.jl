##! Setup
using GLMakie; Makie.inline!(false)
Nx = 100
Ny = 100
xs = range(-0.5, 0.5, Nx)
ys = range(-0.5, 0.5, Ny)
h = step(xs)
dt = 0.5h

##! Setting up figure
fig = Figure(); display(fig)
t = Observable(0.0)
ax = Axis(fig[1, 1], title=@lift("t = $(round($t, sigdigits=3))"))
u = Observable(zeros(Nx, Ny))
heatmap!(ax, u)

##! Defining functions to update state

#* ∂ₜ²u! updates ∂ₜ²u
function ∂ₜ²u!(∂ₜ²u, u, h=h)    # Mutating input 
    for i in 2:(length(xs)-1)
        for j in 2:(length(ys)-1)
            ∂ₜ²u[i-1, j-1] = (u[i+1, j] +   # Takes lines, easy to read.
                            u[i-1, j] + 
                            u[i, j+1] + 
                            u[i, j-1] - 
                           4u[i, j]
            ) / h^2
        end
    end
    return nothing
end

"""
Mutate ∂ₜ²u, ∂ₜu and u, given a stepsize.
"""
function step!(∂ₜ²u, ∂ₜu, u, dt=dt)
    u[2:end-1, 2:end-1] .+= ∂ₜu * dt    #* Update zeroth time derivative using previous slope, like position in kepler
    ∂ₜ²u!(∂ₜ²u, u)                      #* Update ∂ₜ²u, which is same as Δu
    ∂ₜu .= ∂ₜ²u * dt                    #* Update first time derivative using new position, like velocity in kepler
    return nothing
end


##! Init
#* u is Nx × Ny, but as we need points to either side to calculate finite diffs, 
#* the values changing are on the interior.
u_x(x) = cospi(x)
u_y(y) = cospi(y)

# Initial value of u, which will be mutated
u[] = [u_x(x)*u_y(y) for x in xs, y in ys]

∂ₜu = u[][begin+1:end-1, begin+1:end-1]#zeros(size(u[]).-2)
∂ₜ²u = zeros(size(u[]).-2)

#* Testing manual stepping
# display(∂ₜ²u[end÷2, end÷2])
# display(∂ₜu[end÷2, end÷2])
# display(u[][end÷2, end÷2])
# step!(∂ₜ²u, ∂ₜu, u[])

seconds_per_t = 0.1
##! Simulating and plotting
for _ in 1:500
    t0 = time()

    t[] += dt
    step!(∂ₜ²u, ∂ₜu, u[])
    u[] = u[]  # Trigger update
    @show u[][end÷2, end÷2]
    
    sleep(max(0, seconds_per_t*dt - (time()-t0)))
end
##


##! Testing that observables are mutated as expected
function bla!(obs)
    obs[begin+1] = 69
end
blaa = Observable([1, 2, 3])
bla!(blaa[])
blaa[]








using DifferentialEquations
# state u_Δu is 3D array with u in first sheet and Δu in second sheet
u_Δu_0 = []
function my_func(u_Δu_dt, u_Δu, )