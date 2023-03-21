%clear all;
%close all;

Nx = 100;
tmax = 1;
Nt = 120;

dt = tmax / Nt;

xrange = linspace(-0.5, 0.5, Nx);
h = xrange(2) - xrange(1);
v = zeros(1, Nx);
v = exp(-xrange.^2 * 300);
%v = exp(-xrange.^2 * 3000);
%v(40:60) = 1;
%w = zeros(1, Nx-1); 
w = exp(-(xrange(2:end) - h/2).^2 * 300);


for t_idx = 1:Nt
  t = dt * t_idx;

  w = w + dt * diff(v) / h;
  v(2:end-1) = v(2:end-1) + dt / h * diff(w);
  plot(xrange, v);
  axis([min(xrange), max(xrange), -2, 2]);
  usleep(10000);
  drawnow;
end


