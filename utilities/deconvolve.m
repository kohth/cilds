function [DeconvOutput] = deconvolve(RunParam,Observation,varargin)
% Deconvolve takes in fluorescence and uses OASIS to extract spiking
% activity using AR1

% INPUT
%
%      Observation  - Structure containing recorded data (fluorescence
%                     traces)
%                           -Dimensions: 1 x N_TRIAL
%                           -Fields:
%                               y (N_NEURON x T) -- neural data

%% AUTHOR    : Koh Tze Hui
%% DEVELOPED : MATLAB (R2018a)
%% FILENAME  : deconvolve.m
%% LAST CHECKED: 211103 (YYMMDD)

%% =========== Initialize function parameters ==========================
p = inputParser;
p.addRequired('RunParam',@isstruct);
p.addRequired('Data',@isstruct);
p.addParameter('testMask',nan,@islogical);
p.addParameter('trainMask',nan,@islogical);
p.addParameter('splitTestTrain',true,@islogical);
p.parse(RunParam,Observation,varargin{:});


RunParam = p.Results.RunParam;
Observation = p.Results.Data;
trainMask = p.Results.trainMask;
testMask = p.Results.testMask;
splitTestTrain = p.Results.splitTestTrain;

%% === Split testing and training ===
ind = 1:size(Observation,2);
if ~isnan(trainMask)
    RunParam(1).TRAININD = ind(trainMask);
end

if ~isnan(testMask)
    RunParam(1).TESTIND = ind(testMask);
end

if splitTestTrain
    splitTrainInd = RunParam(1).TRAININD;
    splitTestInd = RunParam(1).TESTIND;
    TrainFtrace = Observation(splitTrainInd);
    TestFtrace = Observation(splitTestInd);
else
    TrainFtrace = Observation;
    splitTrainInd = 1:size(Observation,2);
    TestFtrace = [];
end
%% Get parameters
% gamma = convertgamma(Param.GAMMA(gammaIdx),1000/Param.BIN,'toTheirs');
gamma = RunParam(1).GAMMA;
params.p = 1; params.B = 1;

%% === Deconvolve training data ====
for iSplit = 1:size(TrainFtrace,2) % Deconvolution is performed every trial separately
    trainData = num2cell(TrainFtrace(iSplit).y,2); % Format into cell type so can use cellfun
    
    spktrain = zeros(size(TrainFtrace(iSplit).y,1),size(TrainFtrace(iSplit).y,2));
    calcium = zeros(size(TrainFtrace(iSplit).y,1),size(TrainFtrace(iSplit).y,2));
    
    % If Gamma (calcium decay) is given
    if RunParam(1).GAMMAKNOWN
        switch RunParam(1).DECONVOLUTIONTYPE
            case 'threshold'
                decay_time = 0.4; spk_SNR = 0.5; lam_pr = 0.99;
                sMin = cellfun(@(x) spk_SNR*GetSn(x),trainData,'un',0);
                lam = cellfun(@(x) choose_lambda(exp(-1/(1000/RunParam(1).BIN*decay_time)),GetSn(x),lam_pr),trainData,'un',0);
                [cTemp,sTemp,options] = cellfun(@(x,y,z) deconvolveCa(x,...
                    'ar1',gamma,'thresholded','smin',y,'lambda',z,'optimize_b',RunParam(1).OPTIMIZEB),trainData,sMin,lam,'un',0);
            case 'constrained'
                [cTemp,sTemp,options] = cellfun(@(x) deconvolveCa(x,...
                    'ar1',gamma,'constrained','optimize_pars',false,'optimize_b',RunParam(1).OPTIMIZEB),trainData,'un',0);
        end
    else %If Gamma is not given
        switch RunParam(1).DECONVOLUTIONTYPE
            case 'threshold'
                decay_time = 0.4; spk_SNR = 0.5; lam_pr = 0.99;
                sMin = cellfun(@(x) spk_SNR*GetSn(x),trainData,'un',0);
                lam = cellfun(@(x) choose_lambda(exp(-1/(RunParam(1).FRAMERATE*decay_time)),GetSn(x),lam_pr),trainData,'un',0);
                if RunParam(1).INITIALIZEGAMMA
                    [cTemp,sTemp,options] = cellfun(@(x,y,z) deconvolveCa(x,...
                        'ar1',gamma,'thresholded','smin',y,'optimize_pars',RunParam(1).OPTIMIZEPARS,'optimize_b',RunParam(1).OPTIMIZEB,'lambda',z,'thresh_factor',0.99),trainData,sMin,lam,'un',0);
                else
                    [cTemp,sTemp,options] = cellfun(@(x,y,z) deconvolveCa(x,...
                        'ar1','thresholded','smin',y,'optimize_pars',RunParam(1).OPTIMIZEPARS,'optimize_b',RunParam(1).OPTIMIZEB,'lambda',z),trainData,sMin,lam,'un',0);
                end
            case 'constrained'
                if RunParam(1).INITIALIZEGAMMA
                    [cTemp,sTemp,options] = cellfun(@(x) deconvolveCa(x,...
                        'ar1',gamma,'constrained','optimize_pars',RunParam(1).OPTIMIZEPARS,'optimize_b',RunParam(1).OPTIMIZEB),trainData,'un',0);
                    
                else
                    [cTemp,sTemp,options] = cellfun(@(x) deconvolveCa(x,...
                        'ar1','constrained','optimize_pars',RunParam(1).OPTIMIZEPARS,'optimize_b',RunParam(1).OPTIMIZEB),trainData,'un',0);
                end
        end
    end
    gam = zeros(size(options,1),1);
    b = zeros(size(options,1),1);
    smin = zeros(size(options,1),1);
    for i = 1:size(options,1)
        gam(i) = options{i}.pars;
        b(i) = options{i}.b;
        smin(i) = options{i}.smin;
    end
    sTemp = cell2mat(sTemp')';
    cTemp = cell2mat(cTemp')';
    spktrain(:,1:size(sTemp,2)) = sTemp;
    calcium(:,1:size(cTemp,2)) = cTemp;
    
    DeconvOutput(splitTrainInd(iSplit)).b = b;
    DeconvOutput(splitTrainInd(iSplit)).smin = smin;
    DeconvOutput(splitTrainInd(iSplit)).gam = gam;
    DeconvOutput(splitTrainInd(iSplit)).y = spktrain;
    DeconvOutput(splitTrainInd(iSplit)).c = calcium;
end
gam = mean([DeconvOutput(splitTrainInd).gam],2);
b = mean([DeconvOutput(splitTrainInd).b],2);
smin = mean([DeconvOutput(splitTrainInd).smin],2);

%% === Deconvolve testing data using training deconv. parameters ====
for iSplit = 1:size(TestFtrace,2)
    testData = num2cell(TestFtrace(iSplit).y,2);
    spktrain = zeros(size(TrainFtrace(iSplit).y,1),size(TestFtrace(iSplit).y,2));
    calcium = zeros(size(TrainFtrace(iSplit).y,1),size(TestFtrace(iSplit).y,2));
    
    %% Get testing
    switch RunParam(1).DECONVOLUTIONTYPE
        case 'threshold'
            [cTemp,sTemp,options] = cellfun(@(x,y,z,q,p) deconvolveCa(x,...
                'ar1',y,'thresholded','smin',z,'lambda',q,'b',p,'optimize_b',false,'thresh_factor',0.99),testData,num2cell(gam),num2cell(smin),lam,num2cell(b),'un',0);
        case 'constrained'
            [cTemp,sTemp,options] = cellfun(@(x,y,p) deconvolveCa(x,...
                'ar1',y,'constrained','b',p,'optimize_pars',false,'optimize_b',false),testData,num2cell(gam),num2cell(b),'un',0);
            
    end
    gam = zeros(size(options,1),1);
    b = zeros(size(options,1),1);
    smin = zeros(size(options,1),1);
    for i = 1:size(options,1)
        gam(i) = options{i}.pars;
        b(i) = options{i}.b;
        smin(i) = options{i}.smin;
    end
    sTemp = cell2mat(sTemp')';
    cTemp = cell2mat(cTemp')';
    spktrain(:,1:size(sTemp,2)) = sTemp;
    calcium(:,1:size(cTemp,2)) = cTemp;
    DeconvOutput(splitTestInd(iSplit)).b = b;
    DeconvOutput(splitTestInd(iSplit)).smin = smin;
    DeconvOutput(splitTestInd(iSplit)).gam = gam;
    DeconvOutput(splitTestInd(iSplit)).y = spktrain;
    DeconvOutput(splitTestInd(iSplit)).c = calcium;
end

end