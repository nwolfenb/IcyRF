function eps_eff = mixing_shape(eps_e,eps_i,f,N,orientation,model)
% Calculates the effective permittivity based on the permittivity of the
% environment and the volume fraction, permittivity, shape factor,
% orientation and of the inclusion. The dimensionless parameter determines
% the mixing model represented by the following equation:
%
% Syntax:
% eps_eff = mixing_shape(eps_e,eps_i,f,N,orientation,model)
%
% Inputs:
% eps_e         Permittivity of the environment (scalar)
% eps_i         Permittivity of the inclusion (scalar)
% f             Volume fraction of inclusion (scalar or vector)
% N             Shape factor, [Nx Ny Nz] (vector)
% orientation   "aligned" or "random" (string)
% model         "Maxwell Garnett" or "Polder-van Santen"
%
% Outputs:
% eps_eff   Effective permittivity (vector)
%
% Source: 
% Sihvola, A. H. (1999). Electromagnetic mixing formulas and applications
% (No. 47). Iet.
% Shokr (1998) Field Observations and Model Calculations of Dielectric
% Properties of Arctic Sea Ice in the Microwave C-Band
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
    if strcmp(orientation,'aligned')
        eps_eff = eps_e + f*eps_e*(eps_i-eps_e)./(eps_e+(1-f).*N*(eps_i-eps_e));
    elseif strcmp(orientation,'random')
        sigma_num = (eps_i-eps_e)./(eps_e+N*(eps_i-eps_e));
        sigma_den = N*(eps_i-eps_e)./(eps_e+N*(eps_i-eps_e));
        eps_eff = eps_e + eps_e*((f/3)*sum(sigma_num)./(1-(f/3)*sum(sigma_den)));
    else
        error('Input variable orientation must be either "random" or "aligned".')
    end
elseif strcmp(model,'Polder-van Santen')
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
        error('Random Polder-van Santen for randomly orientated inclusions not yet implemented.')
        % Sihvola (1999)
    end
else
    error('Model must be either "Maxwell Garnett" or "Polder-van Santen".')
end
end
