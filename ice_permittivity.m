function eps_ice = ice_permittivity(T,f,sigma)
% Model for the relative permittivity of ice across radio frequencies.
% Combines the models of Gough (1972), Kawada (1978), and Matzler et al.
% (2006).
%
% Syntax:
% eps_ice = ice_permittivity(T,f,sigma)
%
% Inputs:
% T         Temperature (C), scalar or vector
% f         Frequency (Hz), scalar   
% sigma     DC Electrical Conductivity (S/m), scalar or vector
%
% Outputs:
% eps_ice   Relative Permittivity, scalar or vector
%
% Source:
% Matzler et al. (2006)
% Ulaby and Long (2013)
% Kawada (1978)
%
% Author:
% Natalie Wolfenbarger
% nswolfen@gmail.com

%% Convert Temperature to Column Vector
if isrow(T)
    T = T.';
end

%% Celcius to Kelvin
T = T+273.15; % C to K

%% Real Part
eps_real = real(ice_debye(T-273.15,f,sigma));

%% Debye Relaxation
alpha = -imag(ice_debye(T-273.15,f,sigma)).*f/1e9;

%% Infrared Absorption Tail
f = f/1e9; % Hz to GHz
B1 = 0.0207; % K/GHz
b = 335; % K
B2 = 1.16e-11; % GHz^-3
dbeta = exp(-9.963+0.0372*(T-273.16));
beta = (B1./T).*exp(b./T)./(exp(b./T)-1).^2+B2*f.^2+dbeta; % GHz^-1

%% Imaginary Part
eps_imag = alpha./f+beta.*f;

%% Complex Permittivity
eps_ice = eps_real - 1j*eps_imag;

end