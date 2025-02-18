function [alpha,Na] = EMscattering(r,f,epsp,epsb,phi)
% Calculates the attenuation constant and one-way attenutaion rate
% (scattering and absorption) in dB/m due to the presence of spherical
% inclusions.
%
% Syntax:
% [alpha, Na] = EMscattering(r,f,epsp,epsb)
%
% Inputs:
% r         Particle radius (m)
% f         Frequency (Hz)
% epsp      Particle relative permittivity
% epsb      Background relative permittivity
% phi       Porosity
%
% Outputs:
% alpha     Attenuation constant (1/m)
% Na        Attenuation rate (dB/m)
%
% Reference:
% Ulaby and Long (2014)
% Chapter 8
%
% Author:
% Natalie Wolfenbarger
% nswolfen@gmail.com
%
%%
[~, ~, Ee, ~] = Mie_scattering(r, f, epsp, epsb); % Extinction efficiency factor (8-5.2)

Qe = pi*r.^2.*Ee; % Extinction cross section (8-51)
V = 4/3*pi*r.^3; % Volume of a single scatterer
N = phi./V; % Number of scatterers

alpha = (N.*Qe)/2; % Extinction constant
Na = 20*log10(exp(1))*alpha; % One-way attenuation rate (dB/m)
end
