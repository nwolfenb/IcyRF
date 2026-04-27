function eps_s = static_permittivity(T,model,N,flag)
% Models for the static permittivity of ice derived from measurements at
% atmospheric pressure.
%
% Syntax:
% eps_s = static_permittivity(T,model)
% eps_s = static_permittivity(T,model,N)
% eps_s = static_permittivity(T,model,[],flag)
% eps_s = static_permittivity(T,model,N,flag)
%
% Inputs:
% T             Temperature (K), scalar or column vector
% model         Model, string
%               'Auty & Cole (1952)'
%               'Worz & Cole (1969)'
%               'Gough & Davidson (1970)'
%               'Johari & Jones (1975)'
%               'Kawada & Niinuma (1977)'
%               'Kawada (1978)'
%               'Johari & Jones (1978) - parallel'
%               'Johari & Jones (1978) - perpendicular'
%               'Johari & Jones (1978)'
%               'Johari & Whalley (1981)'
%               'Sasaki et al. (2016)'
% N             Number of points defining uncertainty distribution, scalar 
%               [Optional]
% flag          Flag to specify whether to repsect or ignore the range of
%               validity for the model specified, logical
%               true
%               false [Default]
%
% Outputs:
% eps_s         Static permittivity of ice, scalar or column vector
%
% Sources:
% Auty & Cole (1952)
% Worz & Cole (1969)
% Gough & Davidson (1970)
% Johari & Jones (1975)
% Kawada & Niinuma (1977)
% Johari & Jones (1978) - parallel
% Johari & Jones (1978) - perpendicular
% Johari & Jones (1978)
% Johari & Whalley (1981)
% Sasaki et al. (2016)
%
% Author:
% Natalie Wolfenbarger
% nswolfen@gmail.com

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
        T(T<233 | T>273.15) = NaN;
    end
    eps_inf = 3.1;
    C = 37500;
    Tc = -146;

elseif strcmp(model,'Worz & Cole (1969)')
    %% Worz & Cole (1969)
    if flag
        T(T<208 | T>273.15) = NaN;
    end
    eps_inf = 3.1;
    C = 20715;
    Tc = 38;

elseif strcmp(model,'Gough & Davidson (1970)')
    %% Gough & Davidson (1970)
    if flag
        T(T<223 | T>273.15) = NaN;
    end
    eps_inf = 3.1;
    C = 22050;
    Tc = 32;

elseif strcmp(model,'Johari & Jones (1975)')
    %% Johari & Jones (1975)
    if flag
        T(T<109 | T>125) = NaN;
    end
    eps_inf = 3;
    C = 24500;
    Tc = 0;

elseif strcmp(model,'Kawada & Niinuma (1977)')
    %% Kawada & Niinuma (1977)
    if flag
        T(T<123 | T>268) = NaN;
    end
    eps_inf = 3.1;
    C = 22520;
    Tc = 46;
    
elseif strcmp(model,'Kawada (1978)')
    %% Kawada (1978)
    % Estimated from their data by Mastuoka et al. (1996)
    if flag
        T(T<123 | T>273.15) = NaN;
    end
    eps_inf = 3.1;
    C = 23700;
    Tc = 15;

elseif strcmp(model,'Kawada (1978) - parallel')
    %% Kawada (1978)
    if flag
        T(T<123 | T>273.15) = NaN;
    end
    eps_inf = 3.1;
    C = 22500;
    Tc = 46;

elseif strcmp(model,'Kawada (1978) - perpendicular')
    %% Kawada (1978)
    if flag
        T(T<123 | T>273.15) = NaN;
    end
    eps_inf = 3.1;
    C = 23700;
    Tc = 0;

elseif strcmp(model,'Johari & Jones (1978) - parallel')
    %% Johari & Jones (1978) - parallel
    if flag
        T(T<200 | T>271) = NaN;
    end
    eps_inf = 3.2;
    C = 23490;
    sigmaC = 200;
    Tc = 1;
    sigmaTc = 2;

    if exist('N','var')
        if ~isempty(N)
            C = C + sigmaC.*randn(1,N);
            Tc = Tc + sigmaTc.*randn(1,N);
        end
    end

elseif strcmp(model,'Johari & Jones (1978) - perpendicular')
    %% Johari & Jones (1978) - perpendicular
    if flag
        T(T<200 | T>271) = NaN;
    end
    eps_inf = 3.2;
    C = 23260;
    sigmaC = 400;
    Tc = 16;
    sigmaTc = 4;

    if exist('N','var')
        if ~isempty(N)
            C = C + sigmaC.*randn(1,N);
            Tc = Tc + sigmaTc.*randn(1,N);
        end
    end
elseif strcmp(model,'Johari & Jones (1978)')
    %% Johari & Jones (1978)
    if flag
        T(T<200 | T>271) = NaN;
    end
    eps_inf = 3.2;
    C = 23400;
    sigmaC = 200;
    Tc = 15;
    sigmaTc = 2;

    if exist('N','var')
        if ~isempty(N)
            C = C + sigmaC.*randn(1,N);
            Tc = Tc + sigmaTc.*randn(1,N);
        end
    end

elseif strcmp(model,'Johari & Whalley (1981)')
    %% Johari & Whalley (1981)
    if flag
        T(T<133 | T>272) = NaN;
    end
    eps_inf = 3.2;
    C = 24620;
    sigmaC = 201;
    Tc = 6.2;
    sigmaTc = 1.7;

    if exist('N','var')
        if ~isempty(N)
            C = C + sigmaC.*randn(1,N);
            Tc = Tc + sigmaTc.*randn(1,N);
        end
    end

else
    error('Unknown model input.')
end
    eps_s = curieweiss(eps_inf,C,Tc,T);

end

