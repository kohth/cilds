function Param = lds_initializefa(Observation,N_NEURON,N_LATENT,InitParam)
% lds_initializefa uses fa to generate the 
% parameters A, b, R, D, P, h_1, G_1
%                z(t) = Dz(t-1) + u(t)
%                y(t) = Az(t) + b + w(t)
%  where w(t) ~ N(0,R), u(t) ~ N(0,P), z(1)~N(h_1,G_1)
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
%% FILENAME  : lds_initializefa.m
%% LAST CHECKED: 211101 (YYMMDD)

minVarFrac = 0.001;
sumy1 = 0;
N_TRIAL = size(Observation,2);
for iTrial = 1:N_TRIAL
    sumy1 = sumy1 + Observation(iTrial).y(:,1);
    y(:,iTrial) = Observation(iTrial).y(:,1);
end
catdata = concatenate(Observation,'y');

%% === Perform FA ===
[estParam, ~] = fastfa(catdata,N_LATENT);
Z_fa = struct;
for iTrial = 1:N_TRIAL
    FA_result = fastfa_estep(Observation(iTrial).y, estParam);
        
    Z_fa(iTrial).zcurr = FA_result.mean(:,2:end);
    Z_fa(iTrial).zprev = FA_result.mean(:,1:end-1);
    Z_fa(iTrial).zall = FA_result.mean;
    z_1(iTrial,:) = FA_result.mean(:,1)';
end

%% === Set parameters ===
if isfield(InitParam,'A')
    Param.A = InitParam.A;
else
    Param.A = estParam.L; % Use FA loading matrix
end
if isfield(InitParam,'b')
    Param.b = InitParam.b;
else
    Param.b = estParam.d; % Use FA constant
end

if isfield(InitParam,'R')
    Param.R = InitParam.R;
else
    Param.R = diag(estParam.Ph); % Use FA fit noise variance
end

if isfield(InitParam,'D')
    Param.D = InitParam.D;
else
    Param.D = eye(N_LATENT).*0.999; % Heuristic
end

if isfield(InitParam,'P')
    Param.P = InitParam.P;
else
    zcurr = concatenate(Z_fa,'zcurr');
    zprev = concatenate(Z_fa,'zprev');
    Param.P = diag(var(zcurr-Param.D*zprev,0,2)); % Use FA latents to est noise
end

if isfield(InitParam,'h_1')
    Param.h_1 = InitParam.h_1;
else
    Param.h_1 = mean(z_1,1)';
end

if isfield(InitParam,'G_1')
    Param.G_1 = InitParam.G_1;
else
    var_1 = diag(diag(cov(z_1)));
    Param.G_1 = var_1 + diag(diag((z_1-repmat(Param.h_1,1,N_TRIAL)')'*(z_1-repmat(Param.h_1,1,N_TRIAL)')))/N_TRIAL;
    
  if any(diag(Param.G_1) == 0)
      var_all = diag(diag(cov(zcurr')));
      varFloor = minVarFrac*var_all;
      Param.G_1 = max(varFloor, Param.G_1);
      fprintf('Warning: LDS Private variance floor used for one or more observed dimensions in FA.\n');
  end   
end
try
    chol(Param.G_1);
catch
    error('Error. G_1 not positive semi-definite.');
end

end