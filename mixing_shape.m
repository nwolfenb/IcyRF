function eps_eff = mixing_shape(eps_e,eps_i,f,N,orientation,model)
% Calculates the effective permittivity based on the permittivity of the
% environment and the volume fraction, permittivity, shape factor,
% orientation and of the inclusion for a specified mixing model.
%
% Syntax:
% eps_eff = mixing_shape(eps_e,eps_i,f,N,orientation,model)
%
% Inputs:
% eps_e         Environment Permittivity, scalar
% eps_i         Inclusion Permittivity, scalar
% f             Inclusion Volume Fraction, scalar or vector
% N             Shape Factor, vector 
%               Spheres:    [1/3 1/3 1/3]
%               Discs:      [1 0 0]
%               Needles:    [0 1/2 1/2]
% orientation   "aligned" or "random", string
% model         "Maxwell Garnett" or "PVD", string
%
% Outputs:
% eps_eff   Effective Permittivity, scalar or vector
%
% Source:
% Sihvola, A. H. (1999). Electromagnetic mixing formulas and applications
% (No. 47). Iet.
% Shokr, M. E. (1998). Field observations and model calculations of
% dielectric properties of Arctic sea ice in the microwave C-band. IEEE
% transactions on Geoscience and Remote Sensing, 36(2), 463-478.
% 
% Author:
% Natalie Wolfenbarger
% nswolfen@gmail.com
%
%% Check Inputs
if isrow(f)
    f = f';
end

if iscolumn(N)
    N = N';
end

%% Mixing
if strcmp(model,'Maxwell Garnett')
    % Sihvola (1999) 
    if strcmp(orientation,'aligned')
        eps_eff = eps_e + f*eps_e*(eps_i-eps_e)./(eps_e+(1-f).*N*(eps_i-eps_e));
    elseif strcmp(orientation,'random')   
        sigma_num = (eps_i-eps_e)./(eps_e+N*(eps_i-eps_e));
        sigma_den = N*(eps_i-eps_e)./(eps_e+N*(eps_i-eps_e));
        eps_eff = eps_e + eps_e*((f/3)*sum(sigma_num)./(1-(f/3)*sum(sigma_den)));
    else
        error('Input variable orientation must be either "random" or "aligned".')
    end
elseif strcmp(model,'PVD')
    % Shokr (1998)
    if strcmp(orientation,'aligned')
        for n = 1:length(N)
            a = (1-N(n))*ones(size(f));
            b = N(n)*eps_i-eps_e+N(n)*eps_e-f*(eps_i-eps_e);
            c = -N(n)*eps_e*eps_i*ones(size(f));
            eps_roots = cellfun(@roots, num2cell([a b c], 2), 'UniformOutput', false);
            eps_roots = cell2mat(eps_roots').';
            eps_eff(:,n) = max(eps_roots,[],2,'ComparisonMethod','real');
        end
    elseif strcmp(orientation,'random')
        N = mean(N);
        a = (1-N)*ones(size(f));
        b = N*eps_i-eps_e+N*eps_e-f*(eps_i-eps_e);
        c = -N*eps_e*eps_i*ones(size(f));
        eps_roots = cellfun(@roots, num2cell([a b c], 2), 'UniformOutput', false);
        eps_roots = cell2mat(eps_roots').';
        eps_eff = max(eps_roots,[],2,'ComparisonMethod','real');
    end
else
    error('Model must be either "Maxwell Garnett" or "PVD".')
end
end
