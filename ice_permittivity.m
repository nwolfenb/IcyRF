function eps_ice = ice_permittivity(T,f,sigma,model,N,flag)
% Model for the relative permittivity of ice across radio frequencies.
% Combines the models of Gough (1972), Kawada (1978), and Matzler et al.
% (2006).
%
% Syntax:
% eps_ice = ice_permittivity(T,f,sigma)
% eps_ice = ice_permittivity(T,f,sigma,model)
% eps_ice = ice_permittivity(T,f,sigma,model,N)
% eps_ice = ice_permittivity(T,f,sigma,model,N,flag)
%
% Inputs:
% T             Temperature (K), scalar or vector
% f             Frequency (Hz), scalar   
% sigma         DC Electrical Conductivity (S/m), scalar or vector
%               [optional]
% model         Model, structure [Optional]
% .eps_s        Static permittivity model, string
%               'Auty & Cole (1952)'
%               'Worz & Cole (1969)'
%               'Gough & Davidson (1970)'
%               'Johari & Jones (1975)'
%               'Kawada & Niinuma (1977)'
%               'Kawada (1978)' [Default]
%               'Johari & Jones (1978) - parallel'
%               'Johari & Jones (1978) - perpendicular'
%               'Johari & Jones (1978)'
%               'Johari & Whalley (1981)'
%               'Sasaki et al. (2016)'
% .eps_inf      High-frequency permittivity model, string
%               'Gough (1972)' [Default]
% .eps_tau      Relaxation time model, string
%               'Auty & Cole (1952)'
%               'Johari & Jones (1975)'
%               'Johari & Jones (1976) - D2O'
%               'Kawada (1978)' [Default]
%               'Johari & Jones (1978) - parallel'
%               'Johari & Jones (1978) - perpendicular'
%               'Johari & Jones (1978)'
%               'Johari & Whalley (1981)'
%               'Sasaki et al. (2016) - Iha'
%               'Sasaki et al. (2016) - Ihb'
%               'Sasaki et al. (2016) - Ihc'
% N             Number of points defining uncertainty distribution if
%               uncertainty is represented in the selected model, scalar
%               [Optional]
% flag          Flag to specify whether to repsect or ignore the range of
%               validity for the model specified, logical [Optional]
%               true
%               false [Default]
%
% Outputs:
% eps_ice   Relative Permittivity, scalar or vector
%
% Source:
% Matzler et al. (2006)
% Ulaby and Long (2013)
% Kawada (1978)
%
% Author:
% Natalie Wolfenbarger
% nswolfen@gmail.com

%% Check Optional Inputs
% sigma
if ~exist('sigma','var')
    sigma = 0;
end

% model
if ~exist('model','var')
    model.eps_s = 'Kawada (1978)';
    model.eps_inf = 'Gough (1972)';
    model.tau = 'Kawada (1978)';
else
    if isempty(model)
        model.eps_s = 'Kawada (1978)';
        model.eps_inf = 'Gough (1972)';
        model.tau = 'Kawada (1978)';
    end
end

% N
if ~exist('N','var')
    N = [];
end

% flag
if ~exist('flag','var')
    flag = 0;
end

%% Convert Temperature to Column Vector
if isrow(T)
    T = T.';
end

%% Convert Conductivity to Column Vector
if isrow(sigma)
    sigma = sigma.';
end

%% Single Relaxation Debye Model
eps_ice = ice_debye(T,f,sigma,model,N,flag);

%% Real Part
eps_real = real(eps_ice);

%% Debye Relaxation Term (Matzler form)
alpha = -imag(eps_ice).*f/1e9;

%% Infrared Absorption Tail
f = f/1e9; % Hz to GHz
B1 = 0.0207; % K/GHz
b = 335; % K
B2 = 1.16e-11; % GHz^-3
dbeta = exp(-9.963+0.0372*(T-273.16));
beta = (B1./T).*exp(b./T)./(exp(b./T)-1).^2+B2*f.^2+dbeta; % GHz^-1

%% Imaginary Part
eps_imag = alpha./f+beta.*f;

%% Complex Permittivity
eps_ice = eps_real - 1j*eps_imag;

end