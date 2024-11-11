function [alpha, Na] = EMalpha(eps_r,f)
% Calculates the attenuation constant and one-way attenutaion rate in dB/m
% using the relative permittivity.
%
% Inputs:
% eps_r     Relative Permittivity, scalar or vector
% f         Frequency (Hz), scalar or vector
%
% Outputs:
% alpha     Attenuation Constant, scalar or vector
% Na        Attenuation Rate (dB/m), scalar or vector
%
% Reference:
% Ulaby and Long (2014)
% Chapter 2
%
% Author:
% Natalie Wolfenbarger
% nswolfen@gmail.com
%
%% Constants
mu0 = 1.25663706212e-6; % N/A^2
eps0 = 8.8541878128e-12; % F/m

%% Attenuation
w = 2*pi*f;
eps_p = real(eps_r);
eps_pp = -imag(eps_r);
tand = eps_pp./eps_p;

alpha = (w.^2/2*mu0*eps0.*eps_p.*((1+tand.^2).^(1/2)-1)).^(1/2);
Na = 20*log10(exp(1))*alpha;
end