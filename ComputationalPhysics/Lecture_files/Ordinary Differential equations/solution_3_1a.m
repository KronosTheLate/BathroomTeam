clear all;
close all;


N = 180;
t = linspace(0, 100, N);
y = zeros(2, N);
y(:, 1) = [1, 0];
dt = t(2) - t(1);

gamma = 0;	% undamped case
%gamma = .1;	% damped oscillation
%gamma = 2;	% critically damped case
%gamma = 4;	% over-damped case

A = -[0, -1; 1, gamma];

for ii = 2:N
  y(:,ii) = y(:,ii-1) + dt * A * y(:,ii-1);
end

hold on;
plot(t, real(y(1,:)), 'b', 'linewidth', 2)

