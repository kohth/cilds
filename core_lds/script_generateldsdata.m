%% script_generateldsdata
% Generates data from the LDS model
%
%% AUTHOR    : Koh Tze Hui
%% DEVELOPED : MATLAB (R2018a)
%% FILENAME  : script_generateldsdata.m
%% LAST CHECKED: 211101 (YYMMDD)

rng(5);
N_LATENT = 2; N_NEURON = 30; T = 500; N_TRIAL = 500;

%% Generate parameters
TrueParam = lds_initializerandom(N_LATENT,N_NEURON,...
    'coeffR',1,'coeffP',1,'coeffA',0.1);

%% Generate data
[Observation,Latent] = lds_generatedata(TrueParam,'N_TRIAL',N_TRIAL,'T',T);

save('ldsgenerateddata.mat','Observation','Latent','TrueParam');