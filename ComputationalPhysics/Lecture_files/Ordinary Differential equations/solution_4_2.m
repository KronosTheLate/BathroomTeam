clear all;
close all;

tmax = 10000;
expected_convergence = 4;

convergence = [];


gamma = 0;     % undamped case
%gamma = .1;    % damped oscillation
%gamma = 2;     % critically damped case
%gamma = 4;     % over-damped case

% analytical solution:
% we need to identify the critically damped case and handle it specially,
% because it does not follow the form of the non-degenerate solutions.
if (abs(gamma - 2) < 1e-7)
  y_ref_fun = @(t) [(1 + t) .* exp(-t), -t .* exp(-t)];
else
  lambda1 = (-gamma + sqrt(gamma^2 - 4)) / 2;
  lambda2 = (-gamma - sqrt(gamma^2 - 4)) / 2;
  a = lambda2 / (lambda2 - lambda1);
  y_ref_fun = @(t) [...
    a * exp(lambda1 * t) + (1 - a) * exp(lambda2 * t), ...
    lambda1 * a * exp(lambda1 * t) + lambda2 * (1 - a) * exp(lambda2 * t)
    ];
end

A = -[0, -1; 1, gamma];
% 4th-order Taylor expansion:
U_fun = @(A) eye(2) + A * (eye(2) + A / 2 * (eye(2) + A / 3 * (eye(2) + A / 4)));
% 1st-order Taylor expansion:
%U_fun = @(A) eye(2) + A;

% estimate number of steps:
A_norm = max(max(abs(A * tmax)));
initial_error = A_norm^4;
target_error = 1e-8;
N = ceil(log(initial_error / target_error) / log(16));

%N_range = [N];
N_range = 10:30;
NN = length(N_range);
N_range

for N = N_range
  y0 = [1; 0];
  A_working = tmax * A / 2^N;

  U = U_fun(A_working);
  for ii = 1:N
    U = U * U;
  end

  y = U * y0;
  energy = (abs(y(1)).^2 + abs(y(2)).^2) / 2;

  ref_y = y_ref_fun(tmax);
  ref_energy = (abs(ref_y(1)).^2 + abs(ref_y(2)).^2) / 2;

  % amplitude error:
  err = abs(y(1) - ref_y(1));

  % energy error:
  %err = abs(energy - ref_energy);

  convergence = [convergence; err];
end

dt = tmax ./ 2.^N_range;

%loglog(1./N_range, convergence, 'x');
loglog(dt, convergence, 'x');
hold on;
loglog(dt, convergence(end) * (dt / dt(end)).^expected_convergence, '--');
