function Param = cilds_initializelds(Observation,N_NEURON,N_LATENT,InitParam)
%  Function cilds_initializelds uses deconvolution and lds to generate the 
%  parameters G, A, b, Q, B, R, D, P, mu_1, cov_1, h_2, G_2 for the model:
%                z(t) = Dz(t-1) + u(t)
%                c(t) = Gc(t-1) + Az(t) + b + v(t)
%                y(t) = Bc(t) + w(t)
%  where v(t) ~ N(0,Q), w(t) ~ N(0,R), u(t) ~ N(0,P)
%  c(1) ~ N(mu_1,cov_1), z(1)~N(h_2,G_2)
%
%  INPUT:
%      Observation  - Structure containing recorded data (fluorescence
%                     traces)
%                           -Dimensions: 1 x N_TRIAL
%                           -Fields:
%                               y (N_NEURON x T) -- neural data
%
%      N_NEURON     - Scalar indicating number of neurons (dimensionality
%                     of observed data)
%  
%      N_LATENT     - Scalar indicating desired latent dimensionality
%
%      InitParam    - Structure containing the initialization parameters
%  
%% AUTHOR    : Koh Tze Hui
%% DEVELOPED : MATLAB (R2018a)
%% FILENAME  : cilds_initializelds.m
%% LAST CHECKED: 211030 (YYMMDD)

RunParam.N_LATENT = N_LATENT; RunParam.N_NEURON = N_NEURON;

%% === Deconvolve ===
addpath(genpath('../utilities')); % Folder containing deconvolve function
% Initializations for deconvolution
RunParam.GAMMAKNOWN = false;
RunParam.DECONVOLUTIONTYPE = 'constrained';
RunParam.INITIALIZEGAMMA = true;
RunParam.OPTIMIZEB = true;
RunParam.OPTIMIZEPARS = true;
RunParam.TRAININD = 1:size(Observation,2);
if isfield(InitParam,'G')
    RunParam.GAMMA = mean(diag(InitParam.G));
    disp(RunParam.GAMMA);
else
RunParam.GAMMA = 0.9985;
end
[DeconvOutput] = deconvolve(RunParam,Observation,'splitTestTrain',false);

%% === LDS ===
addpath(genpath('../core_lds')); % Access core_lds folder for LDS 
fileHeader = 'ldsinit';
[EstParam] = lds(DeconvOutput, RunParam,'splitTrainTest',false,'fileHeader',...
    fileHeader,'maxIter',100,'partialSaveIter',100);

%% === Get y(1) for mu_1 and cov_1 ===
y_1 = [];
for iTrial = 1:size(Observation,2)    
    y_1 = [y_1 Observation(iTrial).y(:,1)]; 
end

%% === Set parameters ===
if isfield(InitParam, 'G')
    Param.G = InitParam.G;
else
    Param.G = diag(mean([DeconvOutput(:).gam],2));
end
if isfield(InitParam,'D')
    Param.D = InitParam.D;
else
    Param.D = EstParam.D;
end

if isfield(InitParam,'P')
    Param.P = InitParam.P;
else
    Param.P = EstParam.P;
end

if isfield(InitParam,'A')
    Param.A = InitParam.A;
else
    Param.A = EstParam.A;
end

if isfield(InitParam,'b')
    Param.b = InitParam.b;
else
    Param.b = EstParam.b;
end

if isfield(InitParam,'R')
    Param.R = InitParam.R;
else
    c = [DeconvOutput(:).c];
    y = [Observation(:).y];
    Param.R = diag(diag(cov((y-c)')));
end

if isfield(InitParam,'Q')
    Param.Q = InitParam.Q;
else
    Param.Q = EstParam.R;
end


if isfield(InitParam,'B')
    Param.B = InitParam.B;
else
    Param.B = eye(N_NEURON);
end

if isfield(InitParam,'mu_1')
    Param.mu_1 = InitParam.mu_1;
else
    Param.mu_1 = mean(y_1,2);
end

if isfield(InitParam,'cov_1')
    Param.cov_1 = InitParam.cov_1;
else
    var_1 = diag(diag(cov(y_1')));
    cov_1 = var_1 + diag(diag((y_1'-repmat(Param.mu_1,1,size(Observation,2) )')'*(y_1'-repmat(Param.mu_1,1,size(Observation,2) )')))/size(Observation,2) ;
    Param.cov_1 = cov_1;
end
try
    chol(Param.cov_1);
catch
    error('Error. Cov_1 not positive semi-definite.');
end

if isfield(InitParam,'h_1')
    Param.h_2 = InitParam.h_1;
else
    Param.h_2 = EstParam.h_1;
end

if isfield(InitParam,'G_1')
    Param.G_2 = InitParam.G_1;
else
    Param.G_2 = EstParam.G_1;
end
try
    chol(Param.G_2);
catch
    error('Error. G_2 not positive semi-definite.');
end

end