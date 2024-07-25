function eps_ice = ice_gough(T)
% Empirical model for the high frequency relative permittivity of ice as a
% function of temperature.
%
% Syntax:
% eps_ice = ice_permittivity(T)
%
% Inputs:
% T     Temperature (C), scalar or vector
%
% Outputs:
% eps_ice   Relative permittivity 
%
% Source:
% Gough (1972)
% also referenced in Glen and Paren (1975) and Di Paolo et al. (2016)
%
% Author:
% Natalie Wolfenbarger
% nswolfen@gmail.com
%
%%
T = T + 273.15; % K
eps_ice  = 3.093+0.72e-4*T+0.11e-5*T.^2;
end