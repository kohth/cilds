%% Function IHPOISSPROCESS
%  Generates spike times according to an inhomogeneous poisson process using
%  method from Neural Signal Processing notes (Byron Yu)
%
% Input:
%       Param  - Structure containing commonly used parameters
%       lambda - Vector of time varying firing rates
%                   dim(1x1000)
%
% Output:
%       spkTrn - spike train generated from an inhomogeneous poisson
%                process
%                   dim(1x1000)
%
%% AUTHOR    : Koh Tze Hui
%% DEVELOPED : MATLAB (R2017a)
%% FILENAME  : ihpoissprocess.m

function spkTrain = ihpoissprocess(Param,lambda)
lambdaMax = max(lambda);

T = Param.T_GEN*Param.N_TRIAL/Param.N_SPLIT;


%% Get homogeneous spiketimes
nSpike = poissrnd(lambdaMax*T/1000); % N follows a poisson distribution
spkTime = round(sort(T*rand(1,nSpike))); %nearest ms

%% Bin spiketimes into 1ms bins

tVec = 0:Param.T_STEP:T-Param.T_STEP;
spkTrain = histcounts(spkTime,[tVec T]);

%% Remove spikes from homogeneous spiketimes to get inhomogeneous
U = rand(1,size(spkTrain,2));
kept_prob = lambda./lambdaMax;
U(U>kept_prob) = nan;
rem_ind = isnan(U);
spkTrain(rem_ind) = 0;
end
