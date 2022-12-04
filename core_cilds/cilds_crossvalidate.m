function  cilds_crossvalidate(Observation,RunParam,varargin)
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
%
%      zDimList    - Scalar or vector containing latent dimensions to be
%                    used in dimensionality reduction
%                           - Default: [1:5]
%   
%      maxIter     - Scalar indicating maximum number of iterations for EM
%                           - Default: 500
%
%      InitParam   - Structure containing user-defined initialization
%                    parameters
%                           - Possible fields:
%                               A,B,G,D,Q,R,P,b,mu_1,cov_1,h_2,G_2   
%
%      initType    - String to choose initialization type
%                        1. singleInit
%                           1 random draw for all parameters
%                        2. randInit
%                           50 random draws for all parameters
%                        3. ldsInit
%                           Use deconv-LDS to provide an initialization 
%                           point for some parameters
%                        4. fixedInit
%                           All parameters for starting are specified
%
%      FixParam     - Structure indicating which parameters to hold
%                     constant in M step (i.e. don't update the values
%                     beyond initialization)
%                           - Default: empty
%                           - Possible fields:
%                               A,B,G,D,Q,R,P,b,mu_1,cov_1,h_2,G_2
%      fileHeader   - String containing start of file name for saving
%                           - Default: nan
%
%% AUTHOR    : Koh Tze Hui
%% DEVELOPED : MATLAB (R2018a)
%% FILENAME  : cilds_crossvalidate.m
%% LAST CHECKED: 220322 (YYMMDD)

p = inputParser;
p.addRequired('Observation',@isstruct);
p.addRequired('RunParam',@isstruct);
p.addParameter('numFold',2,@isscalar);
p.addParameter('zDimList',[1:5],@isvector);
p.addParameter('maxIter',500,@isscalar);
p.addParameter('InitParam',struct,@isstruct);
p.addParameter('initType','ldsInit',@ischar);
p.addParameter('FixParam',struct,@isstruct);
p.addParameter('fileHeader',nan,@ischar);
p.parse(Observation, RunParam, varargin{:});

Observation = p.Results.Observation;
RunParam = p.Results.RunParam;
numFold = p.Results.numFold;
zDimList = p.Results.zDimList;
maxIter = p.Results.maxIter;
InitParam = p.Results.InitParam;
initType = p.Results.initType;
FixParam = p.Results.FixParam;
fileHeader = p.Results.fileHeader;


% Randomly reorder trials
rand('state', 1);
permFile = strcat(fileHeader,'_cilds_permutedindices.mat');
if ~exist(permFile)
    ind = randperm(RunParam.N_TRIAL);
    save(permFile,'ind');
else
    disp('loading permutation indices');
    load(permFile);
end
Observation = Observation(ind);
% Get division indices for folds
foldDiv = floor(linspace(1, RunParam.N_TRIAL+1, numFold+1));
AllRunParam(1:length(zDimList)) = RunParam;
AllInitParam(1:length(zDimList)) = InitParam;
AllFixParam(1:length(zDimList)) = FixParam;

% Make a structure array of run parameters based on number of evaluated dim
for iDim = 1:length(zDimList)
AllRunParam(iDim).TRAININD = [];
AllRunParam(iDim).TESTIND = [];
AllInitParam(iDim).B = diag(ones(RunParam.N_PROJ,1));
end
RunParam = AllRunParam;
FixParam = AllFixParam;
InitParam = AllInitParam;

for iDim =1:length(zDimList) % Change to parfor if evaluating multiple latent dim
    zDim = zDimList(iDim);
    RunParam(iDim).N_LATENT = zDim;
    fprintf('Processing latent dimensionality = %d\n', zDim);
    for cvf = 1:numFold
        cildsFile = strcat(fileHeader,sprintf('_%.03d_%d',zDim,cvf),...
            '_cilds_result');
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
        
        [CildsParam,CildsResult, LL,~,CildsInitParam] = cilds(Observation,...
            RunParam(iDim),'maxIter',maxIter,'FixParam',FixParam(iDim),'InitParam',...
            InitParam(iDim),'leaveOneOut',...
            true,'fileHeader',cildsFile,'splitTrainTest',isSplit,...
            'initType',initType);
            parsave(CildsResult,CildsParam,strcat(cildsFile,'_test.mat'),...
                LL,CildsInitParam);
        
    end
end
end

