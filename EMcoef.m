function [r,R] = EMcoef(eps_r1,eps_r2)
% Calculates the Fresnel power reflection coefficient given the relative
% permittivity of a two-layered, non-magnetic medium.
%
% Syntax:
% [r, R] = EMcoef(eps_r1,eps_r2)
%
% Inputs:
% eps_r1    Relative Permittivity of first layer, scalar or vector
% eps_r2    Relative Permittivity of second layer, scalar or vector
%
% Outputs:
% R         Power Reflection Coefficient (dB), scalar or vector
% r         Power Reflection Coefficient, scalar or vector
%
% Source:
% Peters et al. (2005)
%
% Author:
% Natalie Wolfenbarger
% nswolfen@gmail.com
%
%%

r = abs((sqrt(eps_r1)-sqrt(eps_r2))./(sqrt(eps_r1)+sqrt(eps_r2))).^2;
R = 10*log10(r);

end