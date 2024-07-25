function eps_eff = debye(eps_s,eps_inf,tau,f,sigma)
% Computes the relative permittivity using the single-relaxation Debye
% model.
%
% Syntax:
% eps_eff = debye(eps_s,eps_inf,tau,f,sigma)
%
% Inputs:
% eps_s     Static permittivity   
% eps_inf   High frequency permittivity  
% tau       Relaxation time (s)
% f         Frequency (Hz)
% sigma     Electrical Conductivity (S/m)
%
% Outputs:
% eps_eff   Relative permittivity
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