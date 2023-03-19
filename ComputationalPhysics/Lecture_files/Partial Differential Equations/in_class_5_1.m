clear all;
%close all;

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

[V, D] = eig(A);
d = diag(D);

d

%plot(imag(d), 'x');


mode = 7;
normalization = max(max(abs(V(1:6, mode)))) * 4;
for t = linspace(0, 30, 300)
  state = real(V(:, mode) * exp(d(mode) * t)) / normalization;
  plot([state(1) - 0.5, state(3), state(5) + 0.5], [state(2), state(4), state(6)], 'o-');
  axis([-1, 1, -1, 1]);
  pause(0.1);
end

