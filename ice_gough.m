function eps_ice = ice_gough(T,N,flag)
% Calculates the high frequency relative permittivity of ice as a function
% of temperature using the empirical model of Gough et al. (1972).
%
% Syntax:
% eps_ice = ice_gough(T)
% eps_ice = ice_gough(T,N)
% eps_ice = ice_gough(T,[],flag)
% eps_ice = ice_gough(T,N,flag)
%
% Inputs:
% T             Temperature (K), scalar or vector
% N             Number of points defining uncertainty distribution, scalar 
%               [Optional]
% flag          Flag to specify whether to repsect or ignore the range of
%               validity for the model specified, logical
%               true
%               false [Default]
%
% Outputs:
% eps_ice       Relative Permittivity, scalar or vector
%
% Source:
% Gough (1972)
% also referenced in Glen and Paren (1975) and Di Paolo et al. (2016)
%
% Author:
% Natalie Wolfenbarger
% nswolfen@gmail.com

%% Check Inputs
if ~exist('flag','var')
    flag = false;
end

%% Check if Column Vector
if isrow(T)
    T = T.';
end

%% Range of Validity
if flag
    T(T>270 | T<2) = NaN;
end

%% Coefficients
A = 3.093;
B = 0.72e-4;
C = 0.11e-5;

%% Uncertainty
if exist('N','var')
    if ~isempty(N)
        sigmaA = 0.003;
        A = A + sigmaA*randn(1,N);
        sigmaB = 0.6e-4;
        B = B + sigmaB*randn(1,N);
        sigmaC = 0.02e-5;
        C = C + sigmaC*randn(1,N);
    end
end

%% Permittivity
eps_ice = A+B.*T+C.*T.^2;

end