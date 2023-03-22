using BenchmarkTools
using GLMakie; Makie.inline!(false)

##! Setting up problem
Nx = 100
Ny = 100
xs = range(-0.5, 0.5, Nx)
ys = range(-0.5, 0.5, Ny)
h = step(xs)
dt = 0.5h

##! Setting up figure
# In function so that I can reset the figure
# without scrolling up
function fig_init(include_3d=true)
    global fig = Figure(fontsize=30)
    display(fig)
    global t = Observable(0.0)
    global mask_active = Observable(false)
    rich_mask_active = @lift($mask_active ? rich("Mask on", color=:green) : rich("Mask off", color=:red))
    rich_time = @lift(rich("t = $(round($t, sigdigits=3))\t\t"))
    title_string = @lift(rich($rich_time, $rich_mask_active))
    global ax = Axis(fig[1, 1], title=title_string)
    ax.aspect=1
    global u = Observable(zeros(Nx, Ny))
    global clims = Observable([-0.5, 0.5])
    global hm = heatmap!(ax, xs, ys, u, colorrange=clims)
    if include_3d
        ax3 = Axis3(fig[2, 1])
        surface!(ax3, xs, ys, u, colorrange=clims)
    end
    Colorbar(fig[:, 2], limits=clims)
    return nothing
end
fig_init()
##! Function definitions
function update_Δu!(u1=u1, Δu=Δu, ∂_x²u1=∂_x²u1, ∂_y²u1=∂_y²u1)
    ∂_x²u1 .= diff(diff(u1, dims=1), dims=1)
    ∂_y²u1 .= diff(diff(u1, dims=2), dims=2)
    Δu .= ∂_x²u1[:, 2:end-1] + ∂_y²u1[2:end-1, :]  #? Updating Δu
    return nothing
end

function update_us!(u0=u0, u1=u1, u2=u2)
    update_Δu!()
    #? Updating the u's
    u2[2:Nx-1, 2:Ny-1] = 2u1[2:Nx-1, 2:Ny-1] - u0[2:Nx-1, 2:Ny-1] + (dt/h)^2*Δu
    u0 .= u1
    u1 .= u2
    return nothing
end
# @btime update_us!()  #¤ Quite fast, very few allocations <3
# display(@benchmark update_us!())
##! Setting up mask
mask_border = [(hypot(x, y) ≤ 0.5 ? 1 : 0) for x in xs, y in ys]
# heatmap(mask_border, axis=(aspect=1,))
function apply_mask!(u0=u0, u1=u1, u2=u2; mask=mask_border)
    u0 .*= mask_border
    u1 .*= mask_border
    u2 .*= mask_border
    return nothing
end
##! Initialization
#* Note that u0, u1 and u2 are used for computation
#* u is used for plotting
u0 = zeros(Nx, Ny)      #! u_old
u1 = zeros(Nx, Ny)      #! u_current
u2 = zeros(Nx, Ny)      #! u_new
∂_x²u1 = zeros(Nx-2, Ny)
∂_y²u1 = zeros(Nx, Ny-2)
Δu = zeros(size(u0).-2)

#? Initial conditions
case_initial_cond = 2
if case_initial_cond == 1
    u1 .= [cospi(x)*cospi(y) for x in xs, y in ys]
    u0 .= u1
elseif case_initial_cond == 2
    u1 .= [exp(-50(x^2+y^2)) for x in xs, y in ys]
    u0 .= u1
else
    error("Initial condition case not recognized.")
end
t[] = 0
u[] = u1  # Trigger update

##! Simulating and plotting
fig_init(false)
# heatmap!(mask_border, colormap=:grays, transparency=true, overdraw=true)
seconds_per_t = 5
for counter in 0:10
    t0 = time()
    t[] += dt
    cond0 = false
    cond1 = true
    cond2 = counter % 100 ≥ 50
    cond3 = counter ≥ 160
    cond4 = 250 ≥ counter ≥ 160
    if cond1
        mask_active[] = true
        apply_mask!()
    else
        mask_active[] = false
    end
    update_us!()
    u[] = log10.(abs.(u1) .+ 1e-14)  # Trigger update
    min_old, max_old = clims[]
    min_now, max_now = extrema(u[])

    contract_limits = true
    if contract_limits
        clims[] = [min_now * 1.05, max_now * 1.05]  # * 1.05 makes a limit of 0.99 include tick at 1
    else
        clims[] = [min(min_old, min_now * 1.05), max(max_old, max_now * 1.05)]  # * 1.05 makes a limit of 0.99 include tick at 1
    end
    
    sleep(max(0, seconds_per_t*dt - (time()-t0)))
end


##!=        Old code below      =!##
##!==============================!##


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