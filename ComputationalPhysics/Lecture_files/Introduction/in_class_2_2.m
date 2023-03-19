clear all;
close all;

deflation = 0.2;
%deflation = 1;

A = [1 2; 3 4];
B = [0 1; -1 0];

x = [1; -1];
%b = [0; 0.5];
b = [0; 1];

history = [];
for ii = 1:1000
  %x = inv(A + norm(x) * B) * b;
  x = (1 - deflation) * x + deflation * inv(A + norm(x) * B) * b;
  history = [history, x];
end

plot(history')
%semilogy(abs(diff(history')))
history(:, end)
