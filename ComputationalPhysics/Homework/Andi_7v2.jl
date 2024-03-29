using LinearAlgebra # Maybe unused
using Plots
# plotly() # Interactive plots
gr()


Nx=100
tmax=300



# h=1
xmin=-0.5*100
xmax=0.5*100
h=(xmax-xmin)/Nx
x=LinRange(xmin:h:xmax)

dt=0.9*h
N=Int(round(tmax/dt))
dt=tmax/N
dt/h

t=LinRange(0:dt:tmax)



# pnow=cos.(π*x)
p_now=exp.(-0.01*x.^2)
p_old=copy(p_now)

dtoverhi2=(dt/h)^2

Nframes=100
frame_nr = 1
for ti in t
    p_new=2*p_now[2:end-1]-p_old[2:end-1]+dtoverhi2*diff(diff(p_now))
    p_old=copy(p_now)
    p_now[2:end-1]=p_new
    if frame_nr <= ti/maximum(t)*Nframes
        frame_nr = frame_nr + 1
        p_now_plot=plot(x,p_now, xlims=(minimum(x), maximum(x)), ylims=(-2, 2))
        display(p_now_plot)
        sleep(0.1)
    end
end


