function eps_eff = mixing(eps_e,eps_i,f,v)
% Calculates the effective permittivity based on the permittivity of the
% environment and the volume fraction and permittivity of the inclusion.
% The dimensionless parameter determines the mixing model represented by
% the following equation:
%
% (eps_eff-eps_e)/(eps_eff+2*eps_e+v*(eps_eff-eps_e)) = 
% f*(eps_i-eps_e)/(eps_i+2*eps_e+v*(eps_eff-eps_e))
%
%
% Syntax:
% eps_eff = mixing(eps_e,eps_i,f,v)
%
% Inputs:
% eps_e     Environment Permittivity, scalar
% eps_i     Inclusion Permittivity, scalar
% f         Inclusion Volume Fraction, scalar or vector
% v         Dimensionless parameter which governs the mixing model, scalar
%           v=0: Maxwell Garnett
%           v=2: Polder-van Santen (Bruggeman)
%           v=3: Coherent Potential
%
% Outputs:
% eps_eff   Effective Permittivity, scalar or vector
%
% Source: 
% Sihvola, A. (2013). Homogenization principles and effect of mixing on
% dielectric behavior. Photonics and Nanostructures-Fundamentals and
% Applications, 11(4), 364-373.
% https://doi.org/10.1016/j.photonics.2013.01.004
%
% Author:
% Natalie Wolfenbarger
% nswolfen@gmail.com
%
%% Check Inputs
if isrow(f)
    f = f';
end

eps_eff = zeros(size(f));

%% Mixing Models
if v == 0 % Maxwell Garnet
    eps_eff = (eps_e.*eps_i+2*eps_e.^2+2*f.*eps_i.*eps_e-2*f.*eps_e.^2)./(eps_i+2*eps_e-f.*(eps_i-eps_e));
elseif v == 2 % Polder-van Santen (Bruggeman)
    eps_eff_mat = [2*ones(size(f)) (eps_i-2*eps_e-3*f.*eps_i+3*f.*eps_e) -eps_i.*eps_e.*ones(size(f))];
    eps_roots = cellfun(@roots, num2cell(eps_eff_mat, 2), 'UniformOutput', false);
    eps_roots = cell2mat(eps_roots').';
    eps_eff = max(eps_roots,[],2,'ComparisonMethod','real');
elseif v == 3 % Coherent Potential
    eps_eff_mat = [3*ones(size(f)) (eps_i-4*eps_e-4*f.*eps_i+4*f.*eps_e) ...
        (-eps_i.*eps_e+eps_e.^2+f.*eps_i.*eps_e-f.*eps_e.^2)];
    eps_roots = cellfun(@roots, num2cell(eps_eff_mat, 2), 'UniformOutput', false);
    eps_roots = cell2mat(eps_roots').';
    eps_eff = max(eps_roots,[],2,'ComparisonMethod','real');
else
    error(['Dimensionless parameter v, must be specified as either',...
        ' v = 0 (Maxwell Garnett), v = 2 (Bruggeman), v = 3 (Coherent Potential)'])
end
end