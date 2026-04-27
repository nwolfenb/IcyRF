function eps_ice = ice_debye(T,f,sigma,model,N,flag)
% Calculates the relative permittivity of ice using the single-relaxation
% Debye model. If Debye parameters are not specified by the user, the model
% adopts the Debye parameters of Kawada (1978) in Matsuoka et al. (1996)
% and Gough et al. (1972).
%
% Syntax:
% eps_ice = ice_debye(T,f,sigma)
% eps_ice = ice_debye(T,f,sigma,model)
% eps_ice = ice_debye(T,f,sigma,model,N)
% eps_ice = ice_debye(T,f,sigma,model,N,flag)
%
% Inputs:
% T             Temperature (K), scalar or vector
% f             Frequency (Hz), scalar   
% sigma         DC electrical conductivity (S/m), scalar or vector
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
% eps_ice   Relative permittivity, scalar or vector
%
% Source:
% see models above
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

%% Static Permittivity
eps_s = static_permittivity(T,model.eps_s,N,flag);

%% High Frequency Permittivity
eps_inf = hf_permittivity(T,model.eps_inf,N,flag);

%% Relaxation Time
tau = relaxation_time(T,model.tau,N,flag);

%% Single-Relaxation Debye Model
eps_ice = debye(eps_s,eps_inf,tau,f,sigma);

end