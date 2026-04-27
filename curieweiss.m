function eps_s = curieweiss(eps_inf,C,Tc,T)
% Calculates the static permittivity of a material using the Curie-Weiss
% model.
%
% Syntax:
% eps_s = curieweiss(eps_inf,C,Tc,T)
%
% Inputs:
% eps_inf   High Frequency Relative Permittivity, scalar or vector  
% C         C?
% Tc        Curie Temperature (K), scalar or vector
% T         Temperature (K), scalar or vector
%
% Outputs:
% eps_s     Static Relative Permittivity, scalar or vector   
%
% Source:
% Kawada & Niinuma (1977)
%
% Author:
% Natalie Wolfenbarger
% nswolfen@gmail.com

%% Static Permittivity
eps_s = eps_inf + C./(T-Tc);

end