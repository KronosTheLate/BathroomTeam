clear all;
close all;

N = 1000;
t = linspace(0, 10, N);
y = zeros(1, N);
y(1) = 1;
dt = t(2) - t(1);
%alpha = -1 + 2 * 1i;
alpha = -1;

for ii = 2:N
  y(ii) = y(ii-1) * (1 + dt * alpha);
end

ref_y = exp(alpha * t);

hold on;
plot(t, real(y), 'b', 'linewidth', 2)
plot(t, real(ref_y), 'r--', 'linewidth', 2);

figure();
loglog(t, abs(y - ref_y) ./ abs(ref_y));
%loglog(t, abs(y - ref_y));
