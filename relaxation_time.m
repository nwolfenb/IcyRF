function tau = relaxation_time(T,model,N,flag)
% Models for the relaxation time of ice across radio frequencies
% derived from measurements at atmospheric pressure.
%
% Syntax:
% tau = relaxation_time(T,model)
% tau = relaxation_time(T,model,N)
% tau = relaxation_time(T,model,[],flag)
% tau = relaxation_time(T,model,N,flag)
%
% Inputs:
% T             Temperature (K), scalar or column vector
% model         Model, string
%               'Auty & Cole (1952)'
%               'Johari & Jones (1975)'
%               'Johari & Jones (1976) - D2O'
%               'Kawada (1978)'
%               'Johari & Jones (1978) - parallel'
%               'Johari & Jones (1978) - perpendicular'
%               'Johari & Jones (1978)'
%               'Johari & Whalley (1981)'
%               'Sasaki et al. (2016) - Iha'
%               'Sasaki et al. (2016) - Ihb'
%               'Sasaki et al. (2016) - Ihc'
% N             Number of points defining uncertainty distribution, scalar 
%               [Optional]
% flag          Flag to specify whether to repsect or ignore the range of
%               validity for the model specified, logical
%               true 
%               false [Default]
%
% Outputs:
% tau       Relaxation time (s), scalar or column vector
%
% Sources:
% Auty & Cole (1952)
% Johari & Jones (1975)
% Kawada (1978)
% Johari & Jones (1978) - parallel
% Johari & Jones (1978) - perpendicular
% Johari & Jones (1978)
% Bitelli et al. (2004)
%
% Author:
% Natalie Wolfenbarger
% nswolfen@gmail.com
%
% Notes: 
% Bitelli et al. (2004) is excluded from this function because their
% published value for B is wrong

%% Check Inputs
if ~exist('flag','var')
    flag = false;
end

%% Check if T is Column Vector
if isrow(T)
    T = T';
end

if strcmp(model,'Auty & Cole (1952)')
    %% Auty & Cole (1952)
    if flag
        T(T<=207 | T>=273.15) = NaN;
    end
    A = 5.3e-16;
    B = 13.25;
    R = 1.987e-3; % kcal/mol/K
    tau = A*exp(B./(R*T));

elseif strcmp(model,'Johari & Jones (1975)')
    %% Johari & Jones (1975)
    if flag
        T(T<=109 | T>=125) = NaN;
    end
    A = 5.0e-20;
    B = 11.3;
    R = 1.987e-3; % kcal/mol/K
    tau = A*exp(B./(R*T));

elseif strcmp(model,'Johari & Jones (1976) - D2O')
    % Note: D2O ice
    if flag
        T(T>273.15 | T<=125) = NaN;
    end
    tau = NaN(size(T));

    % Table 2 in Johari & Whalley (1981)
    Eh = 53.2; % kJ/mol
    Ei = 18.8; % kJ/mol
    El = 46.4; % kJ/mol
    tau0h = 2.36e-15; % s
    tau0i = 5.29e-8; % s
    tau0l = 2.04e-17; % s
    R = 8.314e-3; % kJ/mol/K

    % Temperature Regime Boundaries
    Th = (Eh-Ei)/(R*(log(tau0i)-log(tau0h))); % high-intermediate
    ind1 = find(T>=Th);
    ind2 = find(T<Th);

    % High
    A = tau0h; % s
    B = Eh; % kJ/mol
    tau(ind1) = A*exp(B./(R*T(ind1)));

    % Intermediate
    A = tau0i; % s
    B = Ei; % kJ/mol
    taui = A.*exp(B./(R.*T(ind2)));

    % Low
    A = tau0l; % s
    B = El; % kJ/mol
    taul = A.*exp(B./(R.*T(ind2)));

    % Combined crossover
    x = taui./taul;
    tau(ind2) = 2.*taul./(-x+sqrt(x.^2+4));

elseif strcmp(model,'Kawada (1978)')
    %% Kawada (1978)
    if flag
        T(T<123 | T>=268) = NaN;
    end
    A = 5.3e-16;
    B1 = 5.53e4;
    B2 = 2.26e4;
    R = 8.314; % J/mol/K
    Tcrit = 223; % K
    A2 = A*exp(B1/(R*Tcrit))/exp(B2/(R*Tcrit));
    tau = A*exp(B1./(R*T));
    tau(T<=Tcrit) = A2*exp(B2./(R*T(T<=Tcrit)));

elseif strcmp(model,'Johari & Jones (1978) - parallel')
    % Table III
    if flag
        T(T<=200 | T>=271) = NaN;
    end
    A = 2.3*1e-12; % ps K to K
    Asigma = 3.4*1e-12;
    B = 49.4*1e3; % kJ/mol to J/mol
    Bsigma = 2.5*1e3;
    R = 8.314; % J/mol/K

    if exist('N','var')
        if ~isempty(N)
            A = A + Asigma.*randn(1,N);
            B = B + Bsigma.*randn(1,N);
        end
    end

    tau = (A./T).*exp(B./(R*T));

elseif strcmp(model,'Johari & Jones (1978) - perpendicular')
    % Table III
    if flag
        T(T<=200 | T>=271) = NaN;
    end
    A = 0.23*1e-12; % ps K to K
    Asigma = 0.14*1e-12;
    B = 53.9*1e3; % kJ/mol to J/mol
    Bsigma = 0.8*1e3;
    R = 8.314; % J/mol/K

    if exist('N','var') 
        if ~isempty(N)
            A = A + Asigma.*randn(1,N);
            B = B + Bsigma.*randn(1,N);
        end
    end

    tau = (A./T).*exp(B./(R*T));

elseif strcmp(model,'Johari & Jones (1978)')
    % Table III
    if flag
        T(T<=200 | T>=271) = NaN;
    end
    A = 0.93*1e-12; % ps K to K
    Asigma = 0.22*1e-12;
    B = 51.1*1e3; % kJ/mol to J/mol
    Bsigma = 1.6*1e3;
    R = 8.314; % J/mol/K

    if exist('N','var')
        if ~isempty(N)
            A = A + Asigma.*randn(1,N);
            B = B + Bsigma.*randn(1,N);
        end
    end

    tau = (A./T).*exp(B./(R*T));

elseif strcmp(model,'Johari & Whalley (1981)')
    %% Derived from fit to Johari & Whalley (1981) Figure 5
    % Digitized dataset, so no uncertanity reported

    if flag
        T(T<=166 | T>=257) = NaN;
    end
    A = 5.8e-11;
    B = 28500;
    R = 8.314; % J/mol/K

    tau = A.*exp(B./(R*T));

% elseif strcmp(model,'Bittelli et al. (2004)')
%     %% Bittelli et al. (2004)
%     T(T<=243 | T>=273.15) = NaN;
%     A = 5.3e-16;
%     B = 3154.9;
%     R = 8.314; % J/mol/K
%     tau = A*exp(B./(R*T));

elseif strcmp(model,'Sasaki et al. (2016) - Iha')
    %% Sasaki et al. (2016) - Iha
    % Digitized dataset, so no uncertanity reported
    Tc = 246;
    R = 8.314; % J/mol/K
    if flag
        T(T<=123 | T>=263) = NaN;
    end
    tau = NaN(size(T));

    % High Temperature Regime
    ind1 = find(~isnan(T) & T>=Tc);
    A1 = 1.9e-17;
    B1 = 6.2e+04;
    tau(ind1,:) = A1.*exp(B1./(R*T(ind1)));

    % Low Temperature Regime
    ind2 = find(~isnan(T) & T<=Tc);
    A2 = 6.4e-09;
    B2 = 2.20e+04;
    tau(ind2,:) = A2.*exp(B2./(R*T(ind2)));

elseif strcmp(model,'Sasaki et al. (2016) - Ihb')
    %% Sasaki et al. (2016) - Ihb
    % Digitized dataset, so no uncertanity reported

    Tc = 211;
    R = 8.314; % J/mol/K
    if flag
        T(T<=123 | T>=263) = NaN;
    end
    tau = NaN(size(T));

    % High Temperature Regime
    ind1 = find(~isnan(T) & T>=Tc);
    A1 = 2.4e-13;
    B1 = 4.26e+04;
    tau(ind1,:) = A1.*exp(B1./(R*T(ind1)));

    % Low Temperature Regime
    ind2 = find(~isnan(T) & T<=Tc);
    A2 = 1.5e-06;
    B2 = 1.51e+04;
    tau(ind2,:) = A2.*exp(B2./(R*T(ind2)));

elseif strcmp(model,'Sasaki et al. (2016) - Ihc')
    %% Sasaki et al. (2016) - Ihc
    % Digitized dataset, so no uncertanity reported
    if flag
        T(T<=193 | T>=253) = NaN;
    end
    A = 9.6e-16;
    B = 5.35e4;
    R = 8.314; % J/mol/K

    tau = A.*exp(B./(R*T));

else
    error('Unknown model input.')
end
end
