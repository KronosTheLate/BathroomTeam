clear all;
close all;

%A = [1 2 3 4; 1 5 6 7; 1 1 8 9; 1 1 1 0];
A = [3 4 0; -4 3 0; 0 0 1] / 5;
iter = 100;

N = size(A, 1);

% fix point iteration:
y = rand(N, 1);
lA = [];
change = [];
for ii = 1:iter
  y = y / norm(y);
  x = A * y;
  lambda = y' * (A * y) / (y' * y);
  lA = [lA; lambda];
  change = [change, norm(y - x/norm(x))];
  y = x;
end
y' * (A * y) / (y' * y)


% inverse iteration:
B = inv(A);
y = rand(N, 1);
lB = [];
for ii = 1:iter
  y = y / norm(y);
  x = B * y;
  lambda = y' * (A * y) / (y' * y);
  lB = [lB; lambda];
  y = x;
end
y' * (A * y) / (y' * y)


% shift-and-invert:
mu = -14;
%mu = -1.4;
%mu = 1i;
C = inv(A - mu * eye(N));
y = rand(N, 1);
lC = [];
for ii = 1:iter
  y = y / norm(y);
  x = C * y;
  lambda = y' * (A * y) / (y' * y);
  lC = [lC; lambda];
  y = x;
end
y' * (A * y) / (y' * y)

hold on;
plot(change);
plot(lA, 'ko-')
%plot(lB, 'bo-')
%plot(lC, 'ro-')



