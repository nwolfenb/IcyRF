function eps_inf = hf_permittivity(T,model,N,flag)
% Models for the high frequency permittivity of ice derived from
% measurements of H2O ice at atmospheric pressure.
%
% Syntax:
% eps_inf = hf_permittivity(T,model)
% eps_inf = hf_permittivity(T,model,N)
% eps_inf = hf_permittivity(T,model,[],flag)
% eps_inf = hf_permittivity(T,model,N,flag)
%
% Inputs:
% T             Temperature (K), scalar or column vector
% model         Model, string
%               'Auty & Cole (1952)'
%               'Gough (1972)' [Default]
% N             Number of points defining uncertainty distribution, scalar
%               [Optional]
% flag          Flag to specify whether to repsect or ignore the range of
%               validity for the model specified, logical
%               true
%               false [Default]
%
% Outputs:
% eps_hf        High frequency permittivity, scalar or column vector
%
% Sources:
% Gough (1972)
%
% Author:
% Natalie Wolfenbarger
% nswolfen@gmail.com

%% Check Inputs
if ~exist('model','var')
    model = 'Gough (1972)';
end

if ~exist('N','var')
    N = [];
end

if ~exist('flag','var')
    flag = false;
end

%% Gough (1972)
if strcmp(model,'Gough (1972)')
    eps_inf = ice_gough(T,N,flag);

elseif strcmp(model,'Auty & Cole (1952)')
    eps_inf = 3.1*ones(size(T));
    if flag
        eps_inf(T<(-65+273.15)) = NaN;
    end
else
    error('Unknown model input.')
end

end