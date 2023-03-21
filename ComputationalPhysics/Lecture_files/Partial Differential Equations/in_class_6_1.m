clear all;
%close all;

Nx = 100;
tmax = 300;
Nt = 300;

h = 1;
dt = tmax / Nt

xrange = h * ((1:Nx) - Nx/2);
%p_now = exp(-xrange.^2 * 0.01);
%p_now = exp(-xrange.^2 * 0.1);
%p_now = exp(-xrange.^2 * 0.3);
p_now = zeros(Nx, 1);
p_now(40:60) = 1;
p_past = p_now;


frame_nr = 0;
for t_idx = 1:Nt
  t = dt * t_idx;

  p_new = 2 * p_now(2:end-1) - p_past(2:end-1) + (dt/h)^2 * diff(p_now, 2);
  p_past = p_now;
  p_now(2:end-1) = p_new;
  if frame_nr <= t
    frame_nr = frame_nr + 1;
    plot(xrange, p_now);
    axis([min(xrange), max(xrange), -2, 2]);
    pause(0.0001);
    drawnow;
  end
end


