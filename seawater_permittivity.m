function eps = seawater_permittivity(T,S,f)
% Empirical model for the relative permittivity of seawater,
% represented using a double-Debye model. 
%
% Syntax:
% eps = seawater_permittivity(T,f)
%
% Inputs:
% T     Temperature (K), scalar or vector
% S     Salinity (S), scalar or vector
% f     Frequency (Hz), scalar
%
% Outputs:
% eps   Relative permittivity, scalar or vector or matrix
%
% Source:
% Ulaby and Long (2014)
%
% Author:
% Natalie Wolfenbarger
% nswolfen@gmail.com

%% Kelvin to Celcius
T = T - 273.15;

%% Seawater
eps0 = 8.8541878128e-12;

a1  =  0.46606917e-2;
a2  = -0.26087876e-4;
a3  = -0.63926782e-5;
a4  =  0.63000075e+1;
a5  =  0.26242021e-2;
a6  = -0.42984155e-2;
a7  =  0.34414691e-4;
a8  =  0.17667420e-3;
a9  = -0.20491560e-6;
a10 =  0.58366888e+3;
a11 =  0.12684992e+3;
a12 =  0.69227972e-4;
a13 =  0.38957681e-6;
a14 =  0.30742330e+3;
a15 =  0.12634992e+3;
a16 =  0.37245044e+1;
a17 =  0.92609781e-2;
a18 = -0.26093754e-1;

eps_w0 = 87.85306*exp(-0.00456992*T-a1*S-a2*S.^2-a3*S.*T);
eps_w1 = a4*exp(-a5*T-a6*S-a7*S.*T);
tau1 = (a8+a9*S).*exp(a10./(T+a11))/1e9;
tau2 = (a12+a13*S).*exp(a14./(T+a15))/1e9;
eps_inf = a16+a17*T+a18*S;

sigma0 = 2.903602+8.607e-2*T+4.738817e-4*T.^2-2.991e-6*T.^3+4.3041e-9*T.^4;
P = S*(37.5109+5.45216*S+0.014409*S.^2)./(1004.75+182.283*S+S.^2);
alpha0 = (6.9431+3.2841*S-0.099486*S.^2)./(84.85+69.024*S+S.^2);
alpha1 = 49.843-0.2276*S+0.00198*S.^2;
Q = 1+alpha0*(T-15)./(T+alpha1);
sigma = sigma0.*P.*Q;

eps_p = eps_inf+(eps_w0-eps_w1)./(1+(2*pi*f*tau1).^2)+...
    (eps_w1-eps_inf)./(1+(2*pi*f*tau2).^2);
eps_pp = (2*pi*f*tau1.*(eps_w0-eps_w1))./(1+(2*pi*f*tau1).^2)+...
    (2*pi*f*tau2.*(eps_w1-eps_inf))./(1+(2*pi*f*tau2).^2)+...
    sigma./(2*pi*eps0*f);

eps = eps_p+1i*eps_pp;
end
