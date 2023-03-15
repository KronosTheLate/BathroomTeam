using Symbolics
using Symbolics: derivative
using GLMakie; update_theme!(fontsize=25)
Nx = 1+100
Nt = 1+120

xrange = range(-0.5, 0.5, Nx)
h = step(xrange)
@show h
dt = round(0.9h, sigdigits=5)
println("dt/h = $(round(dt/h, sigdigits=3))")

tmax = 

dt_in_units_of_h = dt/h

@variables x t

u_x = cospi(3x)
# u_x = sinpi(2x)     # Look very cool, though unstable
u_t = cospi(3t)
u = u_x * u_t

u_dt = derivative(u, t)
u_dx = derivative(u, x)

# Initialize v as u_dt and w as u_dx, offset
# by half a step size in each dimension to "stagger" grid.
v = [substitute(u_dt, Dict(x=>xval, t=>0)) for xval in xrange] .|> Symbolics.value
w = [substitute(u_dx, Dict(x=>xval-h/2, t=>-dt/2)) for xval in xrange[2:end]] .|> Symbolics.value

# Observables allow updating plot, updating 
# only the observable.
t_ = Observable(0.0)
v_ = Observable(v)

fig = current_figure()
empty!(fig)
ax1 = Axis(fig[1, 1], title=@lift("dt = $(round(dt/h, sigdigits=3))h\nt = $(round($t_, digits=2))"),
ylabel="∂ₜu", xlabel="x")
y_min, y_max = extrema(v_[])

v_analytical = Observable([substitute(u_dt, Dict(x=>xval, t=>0)) for xval in xrange] .|> Symbolics.value)
errors = Observable(zeros(length(xrange)))
y2_min, y2_max = extrema(errors[])

scatterlines!(ax1, xrange, v_analytical, label="Analytical", color=:green)
scatterlines!(ax1, xrange, v_, label="Numerical", color=:black)
axislegend()

ax2, _ = scatterlines(fig[2, 1], xrange, errors, label="Error", color=:red)
axislegend()
# Uncomment and run to open window. If window is open, 
# the rest of the script only modifies it.
current_figure()

secons_per_t = 2   # a factor in "sleep"

## Live plotting
for _ = 1:Nt
    t0 = time()
    t_[] += dt
    
    w_slope = diff(v_[]) / h
    w += dt * w_slope
    v_slope = diff(w) / h 
    v_[][2:end-1] += dt * v_slope

    v_[] = v_[]  # trigger update of observable
    
    v_analytical[] = [substitute(u_dt, Dict(x=>xval, t=>t_[])) for xval in xrange] .|> Symbolics.value
    errors[] = abs.(v_[] - v_analytical[])
    v_min, v_max = extrema(v_[])
    y_min = min(v_min, y_min)
    y_max = max(v_max, y_max)
    ylims!(ax1, y_min, y_max)

    error_min, error_max = extrema(errors[])
    y2_min, y2_max = min(error_min, y2_min), max(error_max, y2_max)
    ylims!(ax2, y2_min, y2_max)

    sleep(max(0, dt*secons_per_t - (time() - t0)))
end


## Write to file
record(current_figure(), "my_video2.mp4", 1:Nt, 
    framerate=(Nt/tmax/secons_per_t)|>round|>Int) do t_idx
    t_[] = dt * t_idx
    
    global w += dt * diff(v_[]) / h
    v_[][2:end-1] += dt / h * diff(w)
    v_[] = v_[]  # trigger update of observable

    global v_min, v_max = extrema(v_[])
    global y_min = min(v_min, y_min)
    global y_max = max(v_max, y_max)

    ylims!(y_min, y_max)
end

record(current_figure(), "my_video2.mp4", 1:Nt, 
    framerate=(Nt/tmax/secons_per_t)|>round) do t_idx
    t_[] = dt * t_idx
    
    global w += dt * diff(v_[]) / h
    v_[][2:end-1] += dt / h * diff(w)
    v_[] = v_[]  # trigger update of observable

    global v_min, v_max = extrema(v_[])
    global y_min = min(v_min, y_min)
    global y_max = max(v_max, y_max)

    ylims!(y_min, y_max)
end

let framerate = 1.0
    if !isa(framerate, Integer)
        @warn "The given framefrate is not a subtype of `Integer`, and will be rounded to the nearest integer. To supress this warning, provide an integer as the framerate."
        framerate = round(Int, framerate)
    end
end