function [r,R] = EMcoef(eps_r1,eps_r2)
% Computes the Fresnel reflection coefficient given the dielectric
% properties of a two layered medium. The equation assumes the medium is
% non-magnetic (relative permeability = 1).
%
% Syntax:
% [r, R] = EMcoef(eps_r1,eps_r2)
%
% Inputs:
% eps_r1    relative permittivity of first layer
% eps_r2    relative permittivity of second layer
%
% Outputs:
% R         Power reflection coefficient (dB)
% r         Power reflection coefficient
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