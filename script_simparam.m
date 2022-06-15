%% script_simparam
% script creates and saves the run parameters
%% AUTHOR    : Koh Tze Hui
%% DEVELOPED : MATLAB (R2018a)
%% FILENAME  : script_simparam.m

paramFolder = "sim_parameter";
if ~exist(paramFolder)
    mkdir(paramFolder);
end

for simIdx = 4
RunParam = createparameter();
switch simIdx
    case 1 %low noise 20 neurons
        RunParam.NOISE = 0.15;
        RunParam.N_PROJ = 20;
        load('M1_faparams_20.mat');
    case 2 %low noise 50 neurons
        RunParam.NOISE = 0.15;
        RunParam.N_PROJ = 50;
        load('M1_faparams_50.mat');
    case 3 %low noise 94 neurons
        RunParam.NOISE = 0.15;
        load('M1_faparams.mat');
    case 4 %medium noise 20 neurons
        load('M1_faparams_20.mat');
        RunParam.N_PROJ = 20;
    case 5 %medium noise 50 neurons
        load('M1_faparams_50.mat');
        RunParam.N_PROJ = 50;
    case 6 %medium noise 94 neurons
        load('M1_faparams.mat');
    case 7 %high noise 20 neurons
        RunParam.NOISE = 15;
        RunParam.N_PROJ = 20;
        load('M1_faparams_20.mat');
    case 8 %high noise 50 neurons
        RunParam.NOISE = 15;
        RunParam.N_PROJ = 50;
        load('M1_faparams_50.mat');z
    case 9 %high noise 94 neurons
        RunParam.NOISE = 15;
        load('M1_faparams.mat');
    case 10 % medium noise 94 neurons offset 1
        load('M1_faparams.mat');
        faParam.d = ones(size(faParams.d));
    case 11 % medium noise 94 neurons offset 10
        load('M1_faparams.mat');
        faParam.d = 10*ones(size(faParams.d));
    case 12 % medium noise 94 neurons offset 100
        load('M1_faparams.mat');
        faParam.d = 100*ones(size(faParams.d));
end
save(sprintf('./%s/sim%.03d_runparam.mat',paramFolder,simIdx),'faParams','RunParam');
end