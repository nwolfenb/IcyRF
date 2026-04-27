function delta = skindepth(eps_r,f)
% Calculates skin depth of a dielectric material as a function of relative
% permittivity and frequency.
%
% Syntax:
% delta = skindepth(eps_r,f)
%
% Inputs:
% eps_r     Relative permittivity
% f         Frequency (Hz)
%
% Outputs:
% delta     Skin depth (m)
%
% Author:
% Natalie Wolfenbarger
% nswolfen@gmail.com
%
% Reference:
% Ulaby and Long (2014)
% Chapter 2

%% Skin Depth
[alpha, ~] = EMalpha(eps_r,f);
delta = 1./alpha;

end
