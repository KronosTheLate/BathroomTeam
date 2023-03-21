clear all;
close all;


N = 1000;
t = linspace(0, 10, N);
r = zeros(2, N);
p = r;
r(:,1) = [1, 0];
p(:,1) = [0, 1.2];
%p(:,1) = [0, 1];
%p(:,1) = [0, 0.3];

dt = t(2) - t(1);

for ii = 2:N
  r(:,ii) = r(:,ii-1) + dt * p(:,ii-1);
  p(:,ii) = p(:,ii-1) - dt * r(:,ii-1) ./ norm(r(:,ii-1))^3;
end

hold on;
plot(r(1,:), r(2,:), 'b.', 'linewidth', 2)

