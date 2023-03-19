%close all;
%clear all;

%N = 1500;
N = 300;
tmax = 50;
dt = 0.01;
frame_dt = 0.1;
alpha = 1.0 / N;

R = [0, 0, 0];
t_jump = 25;

p = zeros(N, 3);
r = zeros(N, 3);

for ii = 1:3
  r(:, ii) = rand(N, 1) - 0.5;
  p(:, ii) = 3 * (rand(N, 1) - 0.5);
end

frame_nr = 0;
for t = 0:dt:tmax
  if frame_nr * frame_dt <= t
    frame_nr = frame_nr + 1;
    plot(r(:, 1), r(:, 2), 'o');
    hold on
    plot(R(1), R(2), 'rx', 'linewidth', 4, 'markersize', 20);
    hold off;
    axis([-2, 4.0, -3.0, 3.0]);
    title(sprintf('time t=%f', t));
    pause(0.01);
    drawnow;
  end
  dr = r - R;
  F_ext = - (dr(:,1).^2 + dr(:,2).^2 + dr(:,3).^2) .* dr;
  p = p + dt * F_ext;
  r = r + dt * p;
end
