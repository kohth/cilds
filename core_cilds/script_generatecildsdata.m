%% script_generatecildsdata
% Generates data from the CILDS model
%
%% AUTHOR    : Koh Tze Hui
%% DEVELOPED : MATLAB (R2018a)
%% FILENAME  : script_generatecildsdata.m
%% LAST CHECKED: 220322 (YYMMDD)

rng(5);
N_LATENT = 2; N_NEURON = 30; T = 300; N_TRIAL = 100;

InitParam.B = diag(ones(N_NEURON,1));
InitParam.d = zeros(N_NEURON,1);
InitParam.k = zeros(N_LATENT,1);
%% Generate parameters
TrueParam = cilds_initializerandom(N_LATENT,N_NEURON,'InitParam',InitParam,...
    'coeffR',1,'coeffQ',1,'coeffP',1,'coeffA',1);

%% Generate data
[Observation,Latent] = cilds_generatedata(TrueParam,'N_TRIAL',N_TRIAL,'T',T);

save('cildsgenerateddata.mat','Observation','Latent','TrueParam');