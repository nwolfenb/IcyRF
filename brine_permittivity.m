function eps_brine = brine_permittivity(T,f)
% Empirical model for the relative permittivity of seawater brine as a
% function of temperature and frequency. Note that data used to obatin this
% model were limited to temperatures between -5 and -25 C and frequencies
% of 7.5, 9. 11.8, 30, and 40 GHz.
%
% Syntax:
% eps_brine = brine_permittivity(T,f)
%
% Inputs:
% T     Temperature (K), scalar or vector
%
% Outputs:
% eps   Relative permittivity 
%
% Source:
% Stogryn et al. (1985)
%
% Author:
% Natalie Wolfenbarger
% nswolfen@gmail.com

%% Convert Temperature to Column Vector
if isrow(T)
    T = T.';
end

%% Kelvin to Celcius
T = T-273.15;

%% Brine Permittivity
eps_s = (939.66-19.068*T)./(10.737-T); % static permittivity
eps_inf = (82.79+8.19*T.^2)./(15.68 + T.^2); % high frequency permittivity
tau = (1/(2*pi))*(0.1099 + 0.13603e-2*T+0.20894e-3*T.^2+0.28167e-5*T.^3)*1e-9; % relaxation time (s)
sigma = brine_conductivity(T+273.15); % electrical conductivity (S/m)

eps_brine = debye(eps_s,eps_inf,tau,f,sigma);

end