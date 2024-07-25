function eps = water_permittivity(T,f)
% Empirical model for the relative permittivity of pure water as a
% function of temperature and frequency.
%
% Syntax:
% eps = water_permittivity(T,f)
%
% Inputs:
% T     Temperature (C), scalar or vector
% f     Frequency (Hz), scalar of vector
%
% Outputs:
% eps   Relative permittivity 
%
% Source:
% Ulaby and Long (2014)
%
% Author:
% Natalie Wolfenbarger
% nswolfen@gmail.com
%
%%
eps_s = 88.045-0.4147*T+6.295e-4*T.^2+1.075e-5*T.^3; % static permittivity
eps_inf = 4.9; % high frequency permittivity
tau = (1/(2*pi))*(1.1109e-10-3.824e-12*T+6.938e-14*T.^2-5.096e-16*T.^3); % relaxation time (s)
sigma = 0; % electrical conductivity (S/m)
eps = debye(eps_s,eps_inf,tau,f,sigma);

end