using LinearAlgebra # Maybe unused
using Plots
# plotly() # Interactive plots
gr()


Nx = 100
tmax = 300
Nt = 300

h = 1
dt = tmax / Nt

xrange = h * ((1:Nx).- Nx/2)
p_now = exp.(-xrange.^2 * 0.01)
# p_now = exp.(-xrange.^2 * 0.1);
# p_now = exp(-xrange.^2 * 0.3);
# p_now = zeros(Nx, 1);
# p_now(40:60) = 1;
p_past = p_now


# function diff2nd(vector) #should be the same as diff(diff())
#     diff2=zeros(size(vector,1)-2)
#     for i in 2:size(vector,1)-1
#         diff2[i-1]=vector[i+1]-2*vector[i]+vector[i-1]
#     end
#     return diff2
# end


frame_nr = 0
for t_idx = 1:Nt
  t = dt * t_idx

  p_new =  2 * p_now[2:end-1] - p_past[2:end-1] + (dt/h)^2 * diff(diff(p_now))
  p_past = p_now
  p_now[2:end-1] = p_new
  if frame_nr <= t
    frame_nr = frame_nr + 1
    p_plot=Plots.plot(xrange, p_now, xlims=(minimum(xrange), maximum(xrange)), ylims=(-2, 2))
    display(p_plot)
    sleep(0.1)
  end
end









