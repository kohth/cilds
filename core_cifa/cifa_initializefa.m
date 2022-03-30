function Param = cifa_initializefa(Observation,D,P,InitParam)
%  Function cifa_initializefa uses deconvolution and FA to generate the 
%  parameters G, A, b, Q, B, R,mu_1, cov_1 for the model:
%                z(t) ~ N(0,I)
%                c(t) = Gc(t-1) + Az(t) + b + v(t)
%                y(t) = Bc(t) + w(t)
%  where v(t) ~ N(0,Q), w(t) ~ N(0,R)
%  c(1) ~ N(mu_1,cov_1),
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
%% FILENAME  : cifa_initializefa.m
%% LAST CHECKED: 220331(YYMMDD)

minVarFrac = 0.01;
RunParam.N_LATENT = P; RunParam.N_PROJ = D;
if isfield(InitParam,'G')
    RunParam.GAMMA = mean(diag(InitParam.G));
    disp(RunParam.GAMMA);
else
RunParam.GAMMA = 0.9985;
end
RunParam.GAMMAKNOWN = false;
RunParam.DECONVOLUTIONTYPE = 'constrained';
RunParam.INITIALIZEGAMMA = true;
RunParam.OPTIMIZEB = true;
RunParam.OPTIMIZEPARS = true;
RunParam.TRAININD = 1:size(Observation,2);
[DeconvOutput] = deconvolve(RunParam,Observation,'splittesttrain',false);
catdata = concatenate(DeconvOutput,'y');
[EstParam, ~] = fastfa(catdata,P);

c_1 = [];
for iTrial = 1:size(Observation,2)    
    c_1 = [c_1 DeconvOutput(iTrial).c(:,1)]; 
end

%% Set parameters
if isfield(InitParam, 'G')
    Param.G = InitParam.G;
else
    Param.G = diag(mean([DeconvOutput(:).gam],2));
end

if isfield(InitParam,'A')
    Param.A = InitParam.A;
else
    Param.A = EstParam.L;
end

if isfield(InitParam,'b')
    Param.b = InitParam.b;
else
    Param.b = EstParam.d;
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
    Param.Q = diag(EstParam.Ph);
end


if isfield(InitParam,'B')
    Param.B = InitParam.B;
else
    Param.B = eye(D);
end

if isfield(InitParam,'mu_1')
    Param.mu_1 = InitParam.mu_1;
else
    Param.mu_1 = mean(c_1,2);
end

if isfield(InitParam,'cov_1')
    Param.cov_1 = InitParam.cov_1;
else
    var_1 = diag(diag(cov(c_1')));
    cov_1 = var_1 + diag(diag((c_1'-repmat(Param.mu_1,1,size(Observation,2) )')'*(c_1'-repmat(Param.mu_1,1,size(Observation,2) )')))/size(Observation,2) ;
    Param.cov_1 = cov_1;
   if any(diag(Param.cov_1) == 0)
      var_all = diag(diag(cov(c')));
      varFloor = minVarFrac*var_all;
      Param.cov_1 = max(varFloor, Param.cov_1);
      fprintf('Warning: LDS Private variance floor used for one or more observed dimensions in FA.\n');
  end     
end
try
    chol(Param.cov_1);
catch
    error('Error. Cov_1 not positive semi-definite.');
end

end