using Images
w = 800
h = 600
greymin = 0
greymax = 255
xs = 1:w
ys = 1:h

##¤ Sinusoidal grating
λ = 20  # wavelength in pixels?
k = 2π/λ
I_sin = sin.(k*xs)

##¤ Save sine image
let I = I_sin
    I_normalized = (I.-minimum(I))./(maximum(I)-minimum(I)) # make range of I from 0 to 1
    I_greyscale = I_normalized .* (greymax-greymin) .+ greymin # "squezzing" in the effective greyscale range
    I_bmp = I_greyscale./255 # back to 0-1 range for .bmp

    pixelvalues = repeat(transpose(I_bmp), h)
    save("/home/dennishb/Uni/Semester/8. Sem/Advanced Physical Optics/Lab excercise 2/GeneratedGratings/sine_lambda=$λ.bmp", pixelvalues)
end

##¤ Binary grating
pitch = 20
duty = 50  # in percent of pitch
cutoff = pitch * duty/100
zero_to_pitch_vector = [x % pitch for x in xs]  # counts 0 to 20 over and over
I_bin = [(x ≤ cutoff ? 1 : -1) for x in zero_to_pitch_vector]
# I_bin = [ifelse(x%pitch) for x in 2π/pitch * xs]
# [x%pitch for x in 2π/pitch * xs]

##¤ Save binary image
let I = I_bin
    I_normalized = (I.-minimum(I))./(maximum(I)-minimum(I)) # make range of I from 0 to 1
    I_greyscale = I_normalized .* (greymax-greymin) .+ greymin # "squezzing" in the effective greyscale range
    I_bmp = I_greyscale./255 # back to 0-1 range for .bmp

    pixelvalues = repeat(transpose(I_bmp), h)
    save("/home/dennishb/Uni/Semester/8. Sem/Advanced Physical Optics/Lab excercise 2/GeneratedGratings/bin_pitch=$(pitch)_duty=$duty.bmp", pixelvalues)
end



##
#= gratinggen.m below:
% % % gratinggen.m by sebo@mci.sdu.dk
% % % This script generates images 1D gratings of given periodicity, amplitude, size, etc.
clc, clear all, close all;

w=800; % pixels - width of the SLM
h=600; % pixels - height of the SLM

greymin=0; % define grayscale value range
greymax=255;

[X,Y]=meshgrid(1:w,1:h); % 2D grid pixel coordinates

%%% sinusoidal grating
Lambda=20; % grating parameters definitions
kx=2*pi./Lambda; 
I=sin(kx*X); 
% 
% %%% binary grating 
% pitchx=20;
% duty=50; % in %
% I=square(2*pi*pitchx^-1*X,duty);

%%% adjusting amplitude range
I=(I-min(I(:)))./(max(I(:))-min(I(:))); % make range of I from 0 to 1
I=I*(greymax-greymin)+greymin; % "squezzing" in the effective greyscale range
I=I./255; % back to 0-1 range for .bmp



imshow(I);
imwrite(I,'grating.png');

% Fourier Transform 
figure;
FI=fft2(I);
imagesc((fftshift(abs(FI))));
=#