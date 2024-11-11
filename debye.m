function eps_eff = debye(eps_s,eps_inf,tau,f,sigma)
% Calculates the relative permittivity using the single-relaxation Debye
% model.
%
% Syntax:
% eps_eff = debye(eps_s,eps_inf,tau,f,sigma)
%
% Inputs:
% eps_s     Static Relative Permittivity, scalar or vector   
% eps_inf   High Frequency Relative Permittivity, scalar or vector  
% tau       Relaxation Time (s), scalar or vector
% f         Frequency (Hz), scalar
% sigma     DC Electrical Conductivity (S/m), scalar or vector
%
% Outputs:
% eps_eff   Relative Permittivity, scalar or vector
%
% Source:
% Ulaby and Long (2014)
%
% Author:
% Natalie Wolfenbarger
% nswolfen@gmail.com
%
%%
eps_0 = 8.854e-12;
eps_eff = eps_inf +((eps_s-eps_inf)./(1+1j*2*pi*tau*f))-...
    1j*sigma./(2*pi*eps_0*f);
end