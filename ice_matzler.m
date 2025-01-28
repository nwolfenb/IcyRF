function eps_ice = ice_matzler(T,f)
% Calculates the relative permittivity of ice as a function of temperature
% using the semi-empirical model of Matzler (2006).
%
% Syntax:
% eps_ice = ice_matzler(T,f)
%
% Inputs:
% T         Temperature (C), scalar or vector
% f         Frequency (Hz), scalar   
%
% Outputs:
% eps_ice   Complex relative permittivity 
%
% Source:
% Matzler (2006)
% Ulaby and Long (2013)
%
% Author:
% Natalie Wolfenbarger
% nswolfen@gmail.com
%
%%

T = T+273.15; % C to K
f = f/1e9; % Hz to GHz

eps_real = 3.1884+9.1e-4*(T-273.15);

theta = 300./T-1;
alpha = (0.00504+0.0062*theta).*exp(-22.1*theta);

B1 = 0.0207; % K/GHz
b = 335; % K
B2 = 1.16e-11; % GHz^-3
dbeta = exp(-9.963+0.0372*(T-273.16));
beta = (B1./T).*exp(b./T)./(exp(b./T)-1).^2+B2*f.^2+dbeta; % GHz^-1

eps_imag = alpha./f+beta.*f;

eps_ice = eps_real - 1j*eps_imag;

end