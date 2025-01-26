function eps_r = colecole(eps_inf,delta_eps,tau,alpha,f,sigma)
% Calculates the relative permittivity using the Cole-Cole model.
%
% Syntax:
% eps_r = colecole(eps_inf,delta_eps,tau,alpha,f,sigma)
%
% Inputs:
% eps_inf       High Frequency Relative Permittivity, scalar
% delta_eps     Dielectric susceptibility, scalar or row vector  
% tau           Relaxation Time (s), scalar or row vector
% alpha         Cole–Cole distribution parameter, scalar or row vector
% f             Frequency (Hz), scalar or column vector
% sigma         DC Electrical Conductivity (S/m), scalar
%
% Outputs:
% eps_r         Relative Permittivity, scalar or vector
%
% Source:
% Grimm 
%
% Author:
% Natalie Wolfenbarger
% nswolfen@gmail.com
%
%% Input Check
if ~isscalar(f)
    if isrow(f)
        error('Frequency must be either a scalar or a column vector.')
    end
end

if ~isrow(delta_eps) || ~isrow(tau) || ~isrow(alpha)
    error('Dielectric susceptibility, relaxation time, and cole-cole distribution parameter must be either a scalar or a row vector.')
end

%% Cole-Cole
eps_0 = 8.854e-12;
w = 2*pi*f;

eps_r = eps_inf + sum(delta_eps./(1+(1j*w.*tau).^(1-alpha)),2,'omitnan') - ...
    1j*sigma./(eps_0*w);
end