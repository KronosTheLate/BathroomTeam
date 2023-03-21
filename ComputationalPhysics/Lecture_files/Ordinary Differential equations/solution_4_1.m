clear all;
close all;

%Tmax = 1;
Tmax = 15;
N_range = [5 10, 20, 50, 100, 200 500, 1000, 2000, 5000, 10000, 20000, 50000];
NN = length(N_range);
convergence = [];

figure();
for N = N_range
  t = linspace(0, Tmax, N);

  r = zeros(2, N);
  p = r;
  r(:,1) = [1, 0];
  p(:,1) = [0, 1];
  ref_energy = norm(p(:,1))^2 / 2 - 1 / norm(r(:,1));

  dt = t(2) - t(1);

  % Euler:
%  for ii = 2:N
%    p(:,ii) = p(:,ii-1) - dt * r(:,ii-1) ./ norm(r(:,ii-1))^3;
%    r(:,ii) = r(:,ii-1) + dt * p(:,ii-1);
%  end

  % leapfrog:
%  for ii = 2:N
%    p(:,ii) = p(:,ii-1) - dt * r(:,ii-1) ./ norm(r(:,ii-1))^3;
%    r(:,ii) = r(:,ii-1) + dt * p(:,ii);
%  end

  % Runge-Kutta:
  initial_state = [r(:, 1); p(:, 1)];
  operator = @(t, state) [state(3:4); -state(1:2) ./ norm(state(1:2))^3];
  [t, rk_result] = rk4(operator, [0, Tmax], initial_state, dt);
  t = t';
  r = rk_result(:, 1:2)';
  p = rk_result(:, 3:4)';

  energy = (abs(p(1, :)).^2 + abs(p(2, :)).^2) / 2 - 1./sqrt(abs(r(1, :)).^2 + abs(r(2, :)).^2);

  ref_r = [ cos(t); sin(t)];
  %err = sqrt(abs(r(1, :) - ref_r(1, :)).^2 + abs(r(2, :) - ref_r(2, :)).^2);
  err = abs(energy - ref_energy);

  figure(1);
  hold on;
  semilogy(t, err);
  xlabel('time [arb. u.]');
  ylabel('error [arb. u.]');

  figure(2);
  hold on;
  plot(r(1, :), r(2, :));
  xlabel('x-coordinate [arb. u.]');
  ylabel('y-coordinate [arb. u.]');

  convergence = [convergence; err(N)];
end

figure(3);
hold on;
loglog(1./N_range, convergence, 'x');
hold on;
% some slopes as guides to the eye:
loglog(1./N_range, (N_range(1)./N_range), '--');
loglog(1./N_range, (N_range(1)./N_range).^2, '--');
loglog(1./N_range, (N_range(1)./N_range).^4, '--');
xlabel('time step size [arb. u.]');
ylabel('error at end of calculation [arb. u.]');
legend(['calculations'; 'first order (guide to the eye)'; 'second order (guide to the eye)'; 'fourth (guide to the eye)'], 'location', 'southeast');

