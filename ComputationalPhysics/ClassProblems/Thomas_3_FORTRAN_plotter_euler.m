clear
clc
close all


data = load('Thomas_3_euler_output.txt');

time = data(1:1001);
FORTRAN_data = data (1006:2006);

y_ref = exp(-1.*time);


hold on
plot(time,y_ref)
plot(time,FORTRAN_data)
legend('Analytical Value','FORTRAN77 Euler estimate')
title('y = exp(-t)')