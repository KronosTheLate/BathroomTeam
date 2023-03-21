close all;
clear all;

% phase propagation matrix:
P = @(n, D, k) [exp(1i .* n .* D .* k), 0; 0, exp(-1i .* n .* D .* k)];

% interface matrix:
I = @(n1, n2) 0.5 * [1+n1/n2, 1-n1/n2; 1-n1/n2, 1+n1/n2];

% 5 dielectric layer pairs
n_array = [1,2,1,2,1,2,1,2,1,2,1];

% 10 dielectric layer pairs
%n_array = [1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1];

% 20 dielectric layer pairs
%n_array = [1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1];

% 10 dielectric layer pairs with central defect
%n_array = [1,2,1,2,1,2,1,2,1,2,1, 1 ,1,2,1,2,1,2,1,2,1,2,1];

% very weakly lossy central layer
%n_array = [1,2,1,2,1, 2+0.1*1i ,1,2,1,2,1];

% weakly lossy central layer
%n_array = [1,2,1,2,1, 2+5*1i ,1,2,1,2,1];

% mildly lossy central layer
%n_array = [1,2,1,2,1, 2+10*1i ,1,2,1,2,1];

% strongly lossy central layer
%n_array = [1,2,1,2,1, 2+15*1i ,1,2,1,2,1];

N = length(n_array);

k0_range = linspace(0, 3, 5000);

% These arrays will contain the transmission/reflection spectra
T = [];
R = [];
for k0 = k0_range
  Tmat = eye(2);
  % The first N-1 layers:
  for ii = 2:N-1
    Tmat = P(n_array(ii), 1, k0) * I(n_array(ii), n_array(ii-1)) * Tmat;
    %Tmat = Tmat * I(n_array(ii), n_array(ii-1)) * P(n_array(ii), 1, k0);
  end
  % and the last interface:
  Tmat = I(n_array(N), n_array(N-1)) * Tmat;
  %Tmat = Tmat * I(n_array(N), n_array(N-1));

  T = [T, abs(Tmat(1,1) - Tmat(1,2).*Tmat(2,1)./Tmat(2,2))^2];
  R = [R, abs(Tmat(1,2)./Tmat(2,2))^2];
end

figure();
hold on;
plot(k0_range, T, 'bx-');
plot(k0_range, R, 'rx-');
xlabel('vacuum wave number k0');
ylabel('reflectance / transmittance');
legend(['transmittance'; 'reflectance']);
figure();
plot(k0_range, 1-R-T, 'r');
legend(['lost energy']);
xlabel('vacuum wave number k0');
ylabel('energy loss');

