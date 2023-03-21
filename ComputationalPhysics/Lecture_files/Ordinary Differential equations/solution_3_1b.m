clear all;
close all;

tmax = 30;
N_range = [50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000];
NN = length(N_range);
convergence = [];

gamma = 0;     % undamped case
%gamma = .3;    % damped oscillation
%gamma = 2;     % critically damped case
%gamma = 4;     % over-damped case

% analytical solution:
% we need to identify the critically damped case and handle it specially,
% because it does not follow the form of the non-degenerate solutions.
if (abs(gamma - 2) < 1e-7)
  y_ref_fun = @(t) [(1 + t) .* exp(-t); -t .* exp(-t)];
else
  lambda1 = (-gamma + sqrt(gamma^2 - 4)) / 2;
  lambda2 = (-gamma - sqrt(gamma^2 - 4)) / 2;
  a = lambda2 / (lambda2 - lambda1);
  y_ref_fun = @(t) [...
    a * exp(lambda1 * t) + (1 - a) * exp(lambda2 * t); ...
    lambda1 * a * exp(lambda1 * t) + lambda2 * (1 - a) * exp(lambda2 * t)
    ];
end



figure();
for N = N_range
  t = linspace(0, tmax, N);
  y = zeros(2, N);
  A = -[0, -1; 1, gamma];
  y(:, 1) = [1, 0];
  dt = t(2) - t(1);

  for ii = 2:N
    y(:, ii) = y(:, ii-1) + dt * A * y(:, ii-1);
  end
  energy = (abs(y(1, :)).^2 + abs(y(2, :)).^2) / 2;

  ref_y = y_ref_fun(t);
  ref_energy = (abs(ref_y(1, :)).^2 + abs(ref_y(2, :)).^2) / 2;

  % amplitude error:
  err = abs(y(1,:) - ref_y(1, :));

  % energy error:
  %err = abs(energy - ref_energy);

  hold on;
  semilogy(t, err);
  convergence = [convergence; err(N)];
end

figure();
loglog(1./N_range, convergence, 'x', 'linewidth', 2);
hold on;
loglog(1./N_range, convergence(end) * N_range(end) ./ N_range, '--');

