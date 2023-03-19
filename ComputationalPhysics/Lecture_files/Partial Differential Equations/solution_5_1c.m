%close all;
%clear all;

N = 300;
tmax = 50;
dt = 0.01;
frame_dt = 0.05;
alpha = 1.0 / N;

R = [0, 0, 0];
t_jump = 25;

p = zeros(N, 3);
r = zeros(N, 3);

for ii = 1:3
  r(:, ii) = rand(N, 1) - 0.5;
  p(:, ii) = 1 * (rand(N, 1) - 0.5);
end

frame_nr = 0;
for t = 0:dt:tmax
  if frame_nr * frame_dt <= t
    frame_nr = frame_nr + 1;
    plot(r(:, 1), r(:, 2), 'o');
    r_eff = [sum(r(:, 1)), sum(r(:, 2))] / N;
    hold on
    plot(R(1), R(2), 'rx', 'linewidth', 4, 'markersize', 20);
    plot(r_eff(1), r_eff(2), 'r+', 'linewidth', 4, 'markersize', 20);
    hold off;
    axis([-2, 4.0, -3.0, 3.0]);
    title(sprintf('time t=%f', t));
    pause(0.01);
    drawnow;
    if t > t_jump
      R = [1.0, 0, 0];
    end
  end
  dr = r - R;
  F_ext = - (dr(:,1).^2 + dr(:,2).^2 + dr(:,3).^2) .* dr;
  F_Coulomb = zeros(N,3);
  for jj = 1:N
    dr = (r(jj, :) - r(:, :));
    Ftmp = alpha * dr ./ (sqrt(dr(:,1).^2 + dr(:,2).^2 + dr(:,3).^2)).^3;
    Ftmp(jj, :) = 0;
    F_Coulomb(jj, :) = sum(Ftmp, 1);
  end
  p = p + dt * (F_ext + F_Coulomb);
  r = r + dt * p;
end
