function [Tb, Tb_z, Tb1, Tb2, Tb3] = brightness(T,z,eps_r,rs,rb,f,Tsky,phi,r,eps_rp)
% Calculates the microwave brightness temperature using a "cloud" model.
% The cloud model is derived from the radiative transfer equation ignoring
% the source term from scattering and all intermediate reflections,
% including reflections from the top interface.
%
% Syntax:
% [Tb, Tb_z, Tb1, Tb2, Tb3] = brightness(T,z,eps_r,rs,rb,f,Tsky,phi,r,eps_rp)
%
% Inputs:
% T         Temperature (C), vector
% z         Depth (m), vector
% eps_r     Relative permittivity, vector
% rs        Surface power reflection coefficient, scalar
% rb        Basal power reflection coefficient, scalar
% f         Frequency (Hz), scalar   
% Tsky      Sky background temperature, scalar
% phi       Particle volume fraction (porosity), scalar or vector
% r         Particle radius (m), scalar or vector
% eps_rp    Particle relative permittivity, scalar or vector
%
% Outputs:
% Tb        Brightness temperature (K), scalar
% Tb_z      Cumulative brightness temperature (K), vector
%           Defined as the brightness temperature at a depth, z, including
%           the effect of surface transmission
% Tb1       Brightness temperature contribution from medium (K), scalar
% Tb2       Brightness temperature contribution from base (K), scalar
% Tb3       Brightness temperature contribution from reflected downwelling
%           (K), scalar
%
% Source:
% Tan et al. (2015)
% *does not include reflected downwelling term
%
% Author:
% Natalie Wolfenbarger
% nswolfen@gmail.com
%
%% Brightness Temperature

% absorption
[alpha, ~] = EMalpha(eps_r,f);
k_absorption = 2*alpha;

% scattering
if ~exist('phi','var') && ~exist('r','var') && ~exist('eps_rp','var')
    k_scattering = 0;
else
    if length(phi) ~= length(r) || length(phi) ~= length(eps_rp)
        error('Inputs phi, r, and eps_rp must be the same length.')
    else
        k_scattering = zeros(length(phi),1);
        for n = 1:length(phi)
            [Es, ~, ~, ~] = Mie_scattering(r(n), f, eps_rp(n), eps_r(n));
            V_scatter = 4/3*pi*r(n).^3;
            A_scatter = pi*r(n).^2;
            k_scattering(n) = phi./V_scatter.*A_scatter.*Es;
        end
    end
end

% extinction
k_extinction = k_absorption+k_scattering;

% brightness temperature contribution from the medium
Tb_1 = (1-rs)*((T.*k_absorption.*exp(-cumtrapz(z,k_absorption))));
Tb_z1 = -flipud(cumtrapz(flipud(z),flipud(Tb_1)));

% brightness temperature contribution from the base
att_2 = exp(cumtrapz(flipud(z),flipud(k_extinction)));
Tb_z2 = (1-rs)*((1-rb)*T(end)*flipud(att_2));

% downwelling brightness temperature contribution
Tb_1d = flipud((1-rs)*(flipud(T.*k_absorption).*exp(cumtrapz(flipud(z),flipud(k_absorption)))));
Tb_z1d = cumtrapz(z,Tb_1d);
Tb_z3 = (1-rs)*(rb*Tb_z1d(end)*flipud(att_2));

% cumulative brightness temperature
Tb_z = Tb_z1 + Tb_z2 + Tb_z3;

% brightness temperature contribution from the medium
Tb1 = Tb_z1(1);

% brightness temperature contribution from the base
Tb2 = Tb_z2(1);

% brightness temperature contribution from the reflected downwelling
Tb3 = Tb_z3(1);

% brightness temperature
Tb = Tb1 + Tb2 + Tb3 + rs*Tsky;

end


