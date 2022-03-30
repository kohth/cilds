%% Function GAUSSPROCESS
%  Simulates a latent x using a gaussian process, with reference to GPFA
%  Byron Yu, Gaussian Processes for Machine learning and mathematical monk
%  (youtube)
%
% Input:
%       Param - Structure containing commonly used parameters
%       tau   - Scalar value specifying the time scale for the gaussian
%               process
%
% Output:
%       z - generated gaussian process in the form of a 1xn array
%
% Parameters:
%       'kernel' - string specifying type of covariance function to be used
%                       1) 'rbf' (default)
%       'noiseVariance' - scalar specifying noise value during generation

%% AUTHOR    : Koh Tze Hui
%% DEVELOPED : MATLAB (R2017a)
%% FILENAME  : gaussprocess.m

function z = gaussprocess(Param,tau,varargin)
rng('shuffle');
%% Parameters
p = inputParser;
p.addRequired('Param',@isstruct);
p.addRequired('tau',@(x) isnumeric(x) && isscalar(x));
p.addParameter('kernel','rbf',@(x) isstring(x)|| ischar(x));
p.addParameter('noiseVariance',10^-9,@(x) isnumeric(x) && isscalar(x));
p.parse(Param,tau,varargin{:});

Param = p.Results.Param;
tau = p.Results.tau;
kernel = p.Results.kernel; % Selected covariance function
sign = p.Results.noiseVariance;
sigf = 1-sign; % Signal variance
delta = @(x,y) x==y; %Kroneker delta

%% Create covariance function
switch kernel
    case 'rbf'
        k = @(x,y,z) sigf*exp(-((x-y).^2)/(2*z^2)) + sign*delta(x,y); %% Squared exponential
    otherwise
        error('User Defined Error: Unknown covariance kernel specified');
end

if tau < 5000
if size(Param.T_VEC,2) <= 15000
    %% Generate covariance matrix with function
    [C_i,C_j] = meshgrid(Param.T_VEC,Param.T_VEC);
    C = k(C_i,C_j,tau);% Covariance matrix dim(T_GENxT_GEN)
    
    %% Create gaussian process with covariance matrix
    z = mvnrnd(zeros(1,Param.T_GEN*Param.N_TRIAL/Param.N_SPLIT),C);
else % approximate gaussian process using gaussian conditioning
    [C_i,C_j] = meshgrid([1:Param.PIECESIZE+Param.HISTORY],[1:Param.PIECESIZE+Param.HISTORY]);
    C = k(C_i,C_j,tau);% Covariance matrix dim(T_GENxT_GEN)
    temp = mvnrnd(zeros(1,Param.HISTORY),C(1:Param.HISTORY,1:Param.HISTORY));
    A = C(1:Param.HISTORY,1:Param.HISTORY);
    B = C(Param.HISTORY+1:Param.PIECESIZE+Param.HISTORY,Param.HISTORY+1:Param.PIECESIZE+Param.HISTORY);
    D = C(1:Param.HISTORY,Param.HISTORY+1:Param.PIECESIZE+Param.HISTORY);
    sig = B-(D'/A)*D;
    sig = (sig+sig')/2;
    z = temp;
    T = Param.T_GEN*Param.N_TRIAL/Param.N_SPLIT - Param.HISTORY;
    for i = 1:T/Param.PIECESIZE
        mu = (D'/A)*(temp');
        temp = [temp(Param.PIECESIZE+1:Param.HISTORY),mvnrnd(mu,sig)];
        z = [z, temp(Param.HISTORY-Param.PIECESIZE+1:Param.HISTORY)];
        
    end
end
else % Approximate gaussian process through interpolation
    sampleInt = tau/1000;
    [C_i,C_j] = meshgrid([1:Param.PIECESIZE+Param.HISTORY],[1:Param.PIECESIZE+Param.HISTORY]);
    C = k(C_i,C_j,1000);% Covariance matrix dim(T_GENxT_GEN)
    temp = mvnrnd(zeros(1,Param.HISTORY),C(1:Param.HISTORY,1:Param.HISTORY));
    A = C(1:Param.HISTORY,1:Param.HISTORY);
    B = C(Param.HISTORY+1:Param.PIECESIZE+Param.HISTORY,Param.HISTORY+1:Param.PIECESIZE+Param.HISTORY);
    D = C(1:Param.HISTORY,Param.HISTORY+1:Param.PIECESIZE+Param.HISTORY);
    sig = B-(D'/A)*D;
    sig = (sig+sig')/2;
    z = temp;
    T = Param.T_GEN/sampleInt*Param.N_TRIAL/Param.N_SPLIT - Param.HISTORY;
    for i = 1:T/Param.PIECESIZE
        mu = (D'/A)*(temp');
        temp = [temp(Param.PIECESIZE+1:Param.HISTORY),mvnrnd(mu,sig)];
        z = [z, temp(Param.HISTORY-Param.PIECESIZE+1:Param.HISTORY)];
    end    
    TSampled = 1:sampleInt:Param.T_GEN*Param.N_TRIAL/Param.N_SPLIT;
    zSampled = z;
    Tintp = 1:Param.T_GEN*Param.N_TRIAL/Param.N_SPLIT;
    z = interp1(TSampled,zSampled,Tintp);
end
end
