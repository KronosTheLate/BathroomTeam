clear all;
close all;

m1 = 16;
m2 = 12;
m3 = 16;

k1 = 220;
k2 = 18;

K = [k1 0; 0 k2];
M1 = eye(2) / m1;
M2 = eye(2) / m2;
M3 = eye(2) / m3;
O = zeros(2, 2);

A = [...
	O, O, O, M1, O, O; ...
	O, O, O, O, M2, O; ...
	O, O, O, O, O, M3; ...
	-K, K, O, O, O, O; ...
	K, -2*K, K, O, O, O; ...
	O, K, -K, O, O, O; ...
	];

tmax = 30;
fun = @(t, u) A * u;

u0 = [0; 0; 0; 0; 0; 0; -0.3; 1; 0; 0; 0; 0]; 
[t, u] = rk4(fun, [0, tmax], u0, 0.01);

hold on;
plot(t, u(:, 7), 'b');
plot(t, u(:, 8), '--b');

plot(t, u(:, 9), 'r');
plot(t, u(:, 10), '--r');

plot(t, u(:, 11), 'k');
plot(t, u(:, 12), '--k');
hold off;

figure();
for ii = 1:5:length(t)
  X = [ u(ii, 7)-2 u(ii, 9) u(ii, 11)+2];
  Y = [ u(ii, 8) u(ii, 10) u(ii, 12)];
  plot(X, Y, 'o-');
  axis([-3, 3, -3 ,3]);
  drawnow;
end 
