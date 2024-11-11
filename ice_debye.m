function eps_ice = ice_debye(T,f,sigma)
% Single relaxation Debye model for the complex relative permittivity of
% ice
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
%
%%
T = T+273.15; % K

% eps_s
eps_inf_ref = 3.1;
C = 23700;
Tc = 15; % K
eps_s = C./(T-Tc)+eps_inf_ref;

% eps_inf
eps_inf = ice_gough(T-273.15);

% tau
Tcrit = 223; % K

tau0(1) = 5.3e-16;
E(1) = 55.3e3;
E(2) = 22.6e3;
R = 8.314; % J/(mol K)
tau0(2) = tau0(1)*exp(E(1)/(R*Tcrit))/exp(E(2)/(R*Tcrit));
tau = tau0(1)*exp(E(1)./(R*T));
tau(T<=Tcrit) = tau0(2)*exp(E(2)./(R*T(T<=Tcrit)));

eps_ice = debye(eps_s,eps_inf,tau,f,sigma);

end