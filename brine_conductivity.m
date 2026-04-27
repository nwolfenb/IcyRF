function sigma = brine_conductivity(T)
% Empirical model for the electrical conductivity of seawater brine as a
% function of temperature.
%
% Syntax:
% sigma = brine_conductivity(T)
%
% Inputs:
% T     Temperature (K), scalar or vector
%
% Outputs:
% sigma Electrical conductivity (S/m)
%
% Source:
% Stogryn et al. (1985)
%
% Author:
% Natalie Wolfenbarger
% nswolfen@gmail.com

%% Convert Temperature to Column Vector
if isrow(T)
    T = T.';
end

%% Kelvin to Celcius
T = T-273.15;

%% Brine Conductivity
sigma = T;
sigma(T>=-22.9)=-T(T>=-22.9).*exp(0.5193+0.8755e-1*T(T>=-22.9));
sigma(T<-22.9)=-T(T<-22.9).*exp(1.0334+0.11*T(T<-22.9));
% sigma(T<-25) = NaN;
end
