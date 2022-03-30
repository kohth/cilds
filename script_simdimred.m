%% script_simdimred
% Perform different dimensionality reduction methods on fluorescence
%% AUTHOR    : Koh Tze Hui
%% DEVELOPED : MATLAB (R2018a)
%% FILENAME  : script_simdimred.m
%% LAST CHECKED: 220331 (YYMMDD)

timerValue = tic;

simIdx = 4; gammaIdx = 1; tauIdx = 6; %simulation regime
folderName = 'sim_datasample';
redmeth = 'cifa'; % Dimensionality reduction method
maxIter = 1501; %EM iterations

% Load data
FileArr=dir(strcat('./',folderName,sprintf('/sim%.03d*.mat',simIdx)));
load(fullfile('./',folderName,'/',FileArr(1).name));
FileArr = FileArr(2:end);
if ~exist(RunParam.RESULTFOLDER)
    mkdir(RunParam.RESULTFOLDER);
end

latentFile = sprintf("sim%.03d_tau*_sampledlatent.mat",simIdx);
latentArr = dir(strcat('./',folderName,'/',latentFile));
fluoFile = sprintf("sim%.03d_tau*_gamma%d_sampled_fluorescence.mat",simIdx,gammaIdx);
fluoArr = dir(strcat('./',folderName,'/',fluoFile));

%% === Start dimensionality reduction ===
RunParam(1:size(tauCombination,1)) = RunParam;
gam = convertgamma(RunParam(1).GAMMA(gammaIdx),1000/RunParam(1).BIN,'toTheirs');
InitParam.G = diag(gam*ones(RunParam(1).N_PROJ,1));

for iData = 1:size(tauCombination,1)
    data = parload(fullfile('./',folderName,'/',fluoArr(iData).name));
    dataField = fieldnames(data);
    resultFile = strcat('./',RunParam(iData).RESULTFOLDER,'/',...
        sprintf('sim%.03d_tau%.02d_gamma%d',simIdx,iData,gammaIdx));
    if ~exist(resultFile,'file')
        switch redmeth
            case "cilds"
                cilds_crossvalidate(data,RunParam(iData),'zDimList',RunParam(iData).N_LATENT,...
                    'numFold',2,'maxIter',maxIter,'fileHeader',resultFile,...
                    'initParam',InitParam);
            case "lds2stage"
                RunParam(iData).GAMMA = gam;
                lds2stage_crossvalidate(data,RunParam(iData),'zDimList',RunParam(iData).N_LATENT,...
                    'numFold',2,'maxIter',maxIter,'fileHeader',resultFile,'InitParam',InitParam);
            case "lds"
                lds_crossvalidate(data,RunParam(iData),'zDimList',RunParam(iData).N_LATENT,...
                    'numFold',2,'maxIter',maxIter,'fileHeader',resultFile);
            case "cifa"
                cifa_crossvalidate(data,RunParam(iData),'zDimList',RunParam(iData).N_LATENT,...
                    'numFold',2,'maxIter',maxIter,'fileHeader',resultFile,...
                    'initParam',InitParam);
        end
    else
        disp("File already exists");
    end

end
