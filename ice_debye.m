function eps_ice = ice_debye(T,f,sigma)
% Calculates the relative permittivity of ice using the single-relaxation
% Debye model and the Debye parameters of Kawada (1978) in Matsuoka et al.
% (1996) and Gough et al. (1972).
%
% Syntax:
% eps_ice = ice_debye(T,f,sigma)
%
% Inputs:
% T         Temperature (C), scalar or vector
% f         Frequency (Hz), scalar   
% sigma     DC electrical conductivity (S/m), scalar or vector
%
% Outputs:
% eps_ice   Relative permittivity, scalar or vector
%
% Source:
% Kawada (1978) in Matsuoka et al. (1996)
% Gough et al. (1972)
%
% Author:
% Natalie Wolfenbarger
% nswolfen@gmail.com

%% Convert Temperature to Column Vector
if isrow(T)
    T = T.';
end

%% Celcius to Kelvin
T = T+273.15; % K

%% Static Permittivity
eps_inf_ref = 3.1;
C = 23700;
Tc = 15; % K
eps_s = C./(T-Tc)+eps_inf_ref;

%% High Frequency Permittivity
eps_inf = ice_gough(T-273.15);

%% Relaxation Time
Tcrit = 223; % K

tau0(1) = 5.3e-16;
E(1) = 55.3e3;
E(2) = 22.6e3;
R = 8.314; % J/(mol K)
tau0(2) = tau0(1)*exp(E(1)/(R*Tcrit))/exp(E(2)/(R*Tcrit));
tau = tau0(1)*exp(E(1)./(R*T));
tau(T<=Tcrit) = tau0(2)*exp(E(2)./(R*T(T<=Tcrit)));

%% Single-Relaxation Debye Model
eps_ice = debye(eps_s,eps_inf,tau,f,sigma);

end