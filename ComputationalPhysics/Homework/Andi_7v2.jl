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
p_old=p_now

dtoverhi2=(dt/h)^2

for ti in t
    p_new=2*p_now[2:end-1]-p_old[2:end-1]+dtoverhi2*diff(diff(p_now))
    p_old=p_now
    p_now[2:end-1]=p_new
    p_now_plot=plot(x,p_now)
    display(p_now_plot)
    sleep(0.02)
end


