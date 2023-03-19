%clear all;
%close all;

Nx = 300;
tmax = 3;
Nt = 950;

dt = tmax / Nt;

xrange = linspace(-0.25, 0.75, Nx);
h = xrange(2) - xrange(1);
v = zeros(1, Nx);
v = exp(-xrange.^2 * 300);
%v(40:60) = 1;
%w = zeros(1, Nx-1); 
w = -exp(-(xrange(2:end) - h/2).^2 * 300);

eps = ones(1, Nx-2);
eps(round(0.8*Nx):end) = 9;


for t_idx = 1:Nt
  t = dt * t_idx;

  w = w + dt / h * diff(v);
  v(2:end-1) = v(2:end-1) + dt / h * diff(w) ./ eps;
  plot(xrange, v);
  axis([min(xrange), max(xrange), -2, 2]);
  usleep(10000);
  drawnow;
end


