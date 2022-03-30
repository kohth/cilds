function  lds2stage_crossvalidate(Observation,RunParam,varargin)
% Cross-validation to split training of dimensionality reduction model (get
% model parameters) and evaluating of dimensionality reduction model (get
% posteriors z)
%
% INPUT
%
%      Observation  - Structure containing recorded data (fluorescence
%                     traces)
%                           -Dimensions: 1 x N_TRIAL
%                           -Fields:
%                               y (N_NEURON x T) -- neural data
%
%      RunParam    -  Structure containing dimension of observations (no.
%                     neurons) and dimension of latent variables, training
%                     indices and testing indices
%                           -Dimensions: 1 x 1
%                           -Fields:
%                               N_LATENT -- desired latent dimension
%                               TRAININD (only if splitting training and testing
%                                        -- trial indices of training data
%                               TESTIND  (only if splitting training and testing
%                                        -- trial indices of testing data
%
% OPTIONAL INPUT:
%
%      numFold     - Scalar indicating number of crossvalidation folds
%                           - Default: 2
%
%      maxIter     - Scalar indicating maximum number of iterations for EM
%                           - Default: 500
%
%      InitParam   - Structure containing user-defined initialization
%                    parameters
%
%      FixParam     - Structure indicating which parameters to hold
%                     constant in M step (i.e. don't update the values
%                     beyond initialization)
%                           - Default: empty
%                           - Possible fields:
%                               A,D,R,P,b,d,k,h_1,G_1
%
%      fileHeader   - String containing start of file name for saving
%                           - Default: nan
%
%% AUTHOR    : Koh Tze Hui
%% DEVELOPED : MATLAB (R2018a)
%% FILENAME  : lds2stage_crossvalidate.m
%% LAST CHECKED: 220331 (YYMMDD)

p = inputParser;
p.addRequired('Observation',@isstruct);
p.addRequired('RunParam',@isstruct);
p.addParameter('numFold',2,@isscalar);
p.addParameter('zDimList',[1:5],@isvector);
p.addParameter('maxIter',500,@isscalar);
p.addParameter('InitParam',struct,@isstruct);
p.addParameter('FixParam',struct,@isstruct);
p.addParameter('fileHeader',@ischar);
p.parse(Observation, RunParam, varargin{:});

Observation = p.Results.Observation;
RunParam = p.Results.RunParam;
numFold = p.Results.numFold;
zDimList = p.Results.zDimList;
InitParam = p.Results.InitParam;
FixParam = p.Results.FixParam;
maxIter = p.Results.maxIter;
fileHeader = p.Results.fileHeader;

% Randomly reorder trials
rand('state', 1);
permFile = strcat(fileHeader,'_lds2stage_permutedindices.mat');
if ~exist(permFile)
    ind = randperm(RunParam.N_TRIAL);
    save(permFile,'ind');
else
    disp('loading permutation indices');
    load(permFile);
end
Observation = Observation(ind);
foldDiv = floor(linspace(1, RunParam.N_TRIAL+1, numFold+1));

% Make a structure array of run parameters based on number of evaluated dim
AllRunParam(1:length(zDimList)) = RunParam;
AllInitParam(1:length(zDimList)) = InitParam;
AllFixParam(1:length(zDimList)) = FixParam;
for iDim = 1:length(zDimList)
AllRunParam(iDim).TRAININD = [];
AllRunParam(iDim).TESTIND = [];
end
RunParam = AllRunParam;
FixParam = AllFixParam;
InitParam = AllInitParam;

for iDim =1:length(zDimList) % Change to parfor if evaluating multiple latent dim
    
    zDim = zDimList(iDim);
    RunParam(iDim).N_LATENT = zDim;
    fprintf('Processing latent dimensionality = %d\n', zDim);
    for cvf = 0:numFold
        lds2stageFile = strcat(fileHeader,sprintf('_%.03d_%d',zDim,cvf),...
            '_lds2stage_result');
        deconvFile = strcat(fileHeader,...
            sprintf('_%.03d_%d',zDim,cvf),'_deconv.mat');
        if cvf == 0
            fprintf('  Training on all data.\n');
            isSplit = false;
        else
            fprintf('  Cross-validation fold %d of %d.\n', cvf, numFold);
            isSplit = true;
        end
        testMask = false(1, RunParam(iDim).N_TRIAL);
        
        if cvf > 0
            % Set cross-validation folds
            testMask(foldDiv(cvf):foldDiv(cvf+1)-1) = true;
        end
        trainMask = ~testMask;
        
        RunParam(iDim).TRAININD = trainMask;
        RunParam(iDim).TESTIND = testMask;
        
        % Deconvolve data using OASIS
        if ~exist(deconvFile,'file')
        [DeconvOutput] = deconvolve(RunParam(iDim),Observation,'trainMask',trainMask,...
            'testMask',testMask,'splittesttrain',true);
        parsave(DeconvOutput,deconvFile);
        else
            data =parload(deconvFile);
            DeconvOutput = data;
        end
        
        % Use output of deconvolution as input into LDS
        [Lds2stageParam,Lds2stageResult, LL,Lds2stageInitParam] = lds(DeconvOutput,...
            RunParam(iDim),'maxIter',maxIter,'FixParam',FixParam,'InitParam',...
            InitParam,'leaveoneout',...
            true,'fileHeader',lds2stageFile,'splitTrainTest',isSplit,...
            'initType','faInit');
        

            parsave(Lds2stageResult,Lds2stageParam,strcat(lds2stageFile,'_test.mat'),...
                LL,Lds2stageInitParam);
        
    end
end
end

