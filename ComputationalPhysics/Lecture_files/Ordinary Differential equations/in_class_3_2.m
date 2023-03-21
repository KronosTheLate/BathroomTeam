clear all;
close all;


N = 10000;
t = linspace(0, 100, N);
y = zeros(1, N);
y(1) = 1;
dt = t(2) - t(1);
alpha = -1;

for ii = 2:N
  y(ii) = y(ii-1) - dt * (sin(y(ii-1)).^2);
end

ref_y = acot(t + cot(y(1)));

hold on;
semilogy(t, real(y), 'b')
semilogy(t, real(ref_y), 'r--');

%figure();
%semilogy(abs(y - ref_y) ./ abs(ref_y));
