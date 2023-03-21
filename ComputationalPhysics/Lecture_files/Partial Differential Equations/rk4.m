function [T, Y] = rk4(odefun, time_interval, initial, time_step)

% parameters:
%
% odefun:
%   function handle to the right hand side of the ODE
%
% time_interval:
%   2-dimensional vector that contains the starting point and the end point of
%   the time-evolution
%
% initial:
%   initial value condition of the ODE, i.e. a vector of the dimension of the
%   problem that contains y(0).
%
% time_step:
%   size of the time step to be taken to get from time_interval(1) to 
%   time_interval(2). The number of time steps is determined from this.
%
% return:
%
% T:
%   a vector containing the time_interval and all intermediate time steps.
%
% Y:
%   matrix that contains the (vectorial) value of the solution of the ODE at
%   the times listed in T
%
% Example:
%   A = [0, 1; -1, 0];
%   fun = @(t, y) - A * y;
%   y0 = [1; 0];
%   [T, Y] = rk4(fun, [0, 10], y0, 0.01);
%   plot(T, Y(:, 1));

T = time_interval(1):time_step:time_interval(2);
N = length(T);
time_step = (time_interval(2) - time_interval(1)) / N;
T = (time_interval(1):(time_step+1e-14):time_interval(2))';
Y = zeros(N, length(initial));

initial_shape = size(initial);

Y(1, :) = initial;
for ii = 2:N
  t = T(ii);
  k_1 = reshape(time_step * odefun(t, initial), initial_shape(1), initial_shape(2));
  k_2 = reshape(time_step * odefun(t + time_step / 2, initial + k_1 / 2), initial_shape(1), initial_shape(2));
  k_3 = reshape(time_step * odefun(t + time_step / 2, initial + k_2 / 2), initial_shape(1), initial_shape(2));
  k_4 = reshape(time_step * odefun(t + time_step, initial + k_3), initial_shape(1), initial_shape(2));
  initial = initial + (k_1 + 2 * k_2 + 2 * k_3 + k_4) / 6;
  Y(ii, :) = initial;
end

end
