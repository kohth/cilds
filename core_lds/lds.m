function [EstParam, Result, testll,InitParam] = lds(Observation,RunParam,varargin)
%
% Extract neural trajectories using LDS model
%                z(t) = Dz(t-1) + u(t)
%                y(t) = Az(t) + b + w(t)
%  where w(t) ~ N(0,R), u(t) ~ N(0,P), z(1)~N(h_1,G_1)
%
% Example use:
%      1. [EstParam, Result] =
%               lds(Observation,RunParam);
%      2. [EstParam, Result, testll,InitParam] =
%               lds(Observation,RunParam,'maxIter',1000);
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
% OUTPUT
%
%      EstParam     - Structure containing estimated model parameters from
%                     maximization step
%                           -Dimensionality: typically 1 x 1, but if parameters are saved
%                            from all iterations, then 1 x maxIter
%                           -Fields:
%                               A (N_NEURONxN_LATENT),
%                               D (N_LATENTxN_LATENT),
%                               R (N_NEURONxN_NEURON),
%                               P (N_NEURONxN_LATENT),
%                               b (N_NEURONx1),
%                               h_1 (N_LATENTx1),
%                               G_1 (N_LATENTxN_LATENT)
%
%      Result       - Structure containing estimated posteriors from
%                     expectation step
%                           -Dimensions: 1 x N_TRIAL
%                           -Fields:
%                               z, c. If leaveoneout toggled, flProj and frProj saved
%                               (predicted fluorescence and firing rate)
%
%      testll       - Vector containing values of log-likelihood for test
%                     data
%                           -Dimensions: 1 x maxIter
%
%      InitParam    - Structure containing the initialization parameters
%
% OPTIONAL INPUT:
%
%      maxIter     - Scalar indicating maximum number of iterations for EM
%                           - Default: 500
%
%      InitParam   - Structure containing user-defined initialization
%                    parameters
%                           -Dimensionality: 1x1
%                           -Possible fields:
%                               A (N_NEURONxN_LATENT),
%                               D (N_LATENTxN_LATENT),
%                               R (N_NEURONxN_NEURON),
%                               P (N_NEURONxN_LATENT),
%                               b (N_NEURONx1),
%                               h_1 (N_LATENTx1),
%                               G_1 (N_LATENTxN_LATENT)
%
%      initType    - String to choose initialization type
%                        1. singleInit
%                           1 random draw for all parameters
%                        2. randInit
%                           50 random draws for all parameters
%                        3. faInit
%                           Use FA to provide an initialization point for
%                           some parameters
%                        4. fixedInit
%                           All parameters for starting are specified
%
%      FixParam     - Structure indicating which parameters to hold
%                     constant in M step (i.e. don't update the values
%                     beyond initialization)
%                           - Default: empty
%                           - Possible fields:
%                               A,D,R,P,b,d,k,h_1,G_1
%
%    splitTrainTest - Boolean indicating whether to split training and
%                     testing set data. If true, use training and testing
%                     indices defined in RunParam
%                           - Default: false
%
%      leaveOneOut  - Boolean indicating whether to predict held out
%                     fluorescence and firing rates
%                           - Default: false
%
%  partialSaveIter  - Scalar indicating how often to save data during
%                     fitting
%                           - Default: 300
%
%      fileHeader   - String containing start of file name for saving
%                           - Default: nan
%
%% AUTHOR    : Koh Tze Hui
%% DEVELOPED : MATLAB (R2018a)
%% FILENAME  : lds.m
%% LAST CHECKED: 211101 (YYMMDD)

%% =========== Initialize function parameters ==========================
p = inputParser;
p.addRequired('Observation',@isstruct);
p.addRequired('RunParam',@isstruct);
p.addParameter('maxIter',500,@isscalar);
p.addParameter('InitParam',struct, @isstruct);
p.addParameter('initType','faInit',@ischar);
p.addParameter('FixParam',struct, @isstruct);
p.addParameter('splitTrainTest',true,@islogical);
p.addParameter('leaveOneOut',false,@islogical);
p.addParameter('partialSaveIter',300,@isscalar);
p.addParameter('fileHeader',nan,@ischar);
p.parse(Observation, RunParam,varargin{:});


Observation = p.Results.Observation;
RunParam = p.Results.RunParam;
maxIter = p.Results.maxIter;
InitParam = p.Results.InitParam;
initType = p.Results.initType;
FixParam = p.Results.FixParam;
splitTrainTest = p.Results.splitTrainTest;
leaveOneOut = p.Results.leaveOneOut;
partialSaveIter = p.Results.partialSaveIter;
fileHeader = p.Results.fileHeader;

%% ========== Split data into training and testing =============
latDim = RunParam.N_LATENT;
obsDim = size(Observation(1).y,1);

if splitTrainTest
    trainInd = RunParam.TRAININD;
    testInd = RunParam.TESTIND;
    TrainObs = Observation(trainInd);
    TestObs = Observation(testInd);
else
    TrainObs = Observation;
    TestObs = TrainObs;
end

%% ======== Initialize model parameters ============
% Create a filename using the date if none provided
if isnan(fileHeader)
    formatOut = 'yymmdd';
    currDate = datestr(now,formatOut);
    fileHeader = strcat(currDate,'_cilds');
end

% Initialize
if ~exist(strcat(fileHeader,'_partial.mat'),'file') % Initialization needed only if starting a new fit
    switch initType
        case 'randInit'
            InitParam = lds_generate_param(latDim,obsDim,'InitParam',InitParam);
            prevll = -inf;
            for iIter = 1:50
                [~,~,~,currll] = lds_estep(InitParam,TrainObs);
                if currll > prevll
                    EstParam = InitParam;
                    prevll = currll;
                end
                InitParam = lds_initializerandom(latDim,obsDim,'InitParam',InitParam);
            end
        case 'fixedInit'
            EstParam = InitParam;
        case 'faInit'
            EstParam = lds_initializefa(TrainObs,obsDim,latDim,InitParam);
        otherwise
            error('Error. Unknown initialization method.');
    end
    InitParam = EstParam;
end


%% ============= Training phase ================
tol = 1e-13; % tolerance point for determining convergence of EM
saveIter = 1;
ll = zeros(maxIter+1,1);  % training loglikelihood, should be monotonically increasing (barring machine error)

% Load existing EM file to continue from where algorithm was left off
isLoaded = false; % we have yet to load the EM file
if exist(strcat(fileHeader,'_partial.mat'),'file') %if it exists, load it
    load(strcat(fileHeader,'_partial.mat'));
    disp('loading partial results');
    startIter = iIter; %which iteration did we leave off at
    if startIter > maxIter
        error('Warning. loaded iterations more than specified for run');
    end
    isLoaded = true; % we have loaded the EM file
else
    % Start of a new EM run with one e-step
    [xf,vf,vj,ll(1)] = lds_estep(EstParam,TrainObs); % CILDS expectation step
    startIter = 2; % EM iteration
end

if partialSaveIter < maxIter % save the parameters
    SaveParam(1) = EstParam;
end

%% Start expectation maximization algorithm
if maxIter>1
    for iIter = startIter:maxIter+1
        dispText = fprintf('LDS EM iteration %d of %d\n',iIter-1,maxIter);
        if ~isLoaded
            % ====== Maximization Step ======
            EstParam = lds_mstep(xf,vf,vj,EstParam,TrainObs,FixParam);
        else
            isLoaded = false; % set the algorithm back to unloaded state so that m step can be performed next iteration
        end
        
        % ====== Expectation Step ======
        [xf,vf,vj,ll(iIter)] = lds_estep(EstParam,TrainObs);
        
        % ===== Stopping criterion ======
        % Stop EM algorithm if loglikelihood decreases
        if ll(iIter)-ll(iIter-1) < 0
            error('Warning. Loglikelihood decreased.');
        end
        
        % Stop EM algorithm if the ll change is within a certain tolerance
        if (ll(iIter)-ll(1))<(1+tol)*(ll(iIter-1)-ll(1))
            break;
        end
        
        % ==== Save the parameters every N iterations, default 300 =====
        if partialSaveIter < maxIter && rem(iIter-1,partialSaveIter)==0 && iIter <= maxIter+1
            SaveParam(iIter) = EstParam;
            saveIter = saveIter+1;
            if ~exist('SaveParam')
                SaveParam = nan;
            end
            save(strcat(fileHeader,'_partial.mat'),'EstParam','ll','iIter','SaveParam','InitParam');
        end
        fprintf(repmat('\b',1,dispText))
    end
end

dispText = fprintf('LDS EM iteration %d of %d\n',iIter-1,maxIter);

%% ============= Testing phase ================
[xf,~,~,testll]  = lds_estep(EstParam,TestObs);
z = xf(:,1:latDim,:); % Extract low-dimensional neural trajectories
Result = struct;
for iTrial = 1:size(z,1)
    Result(iTrial).z = squeeze(z(iTrial,:,:));
end

if leaveOneOut
    % Iteratively leave out one neuron, recompute latent estimates, and
    % reconstruct fluorescence of left-out neuron
    fr_proj = zeros(size(z,1),obsDim,size(z,3));
    fl_proj = zeros(size(z,1),obsDim,size(z,3));
    for iDim = 1:obsDim
        LooObservation = struct;
        mi = [1:(iDim-1) (iDim+1):obsDim];
        for nTrial = 1:size(TestObs,2)
            LooObservation(nTrial).y = TestObs(nTrial).y(mi,:);
        end
        LooEstParam = looparam(EstParam,mi);
        [zLoo,~,~] = lds_estep(LooEstParam,LooObservation);
        
        [fr_proj(:,iDim,:),fl_proj(:,iDim,:)] =loof(iDim,EstParam,zLoo,size(TestObs,2),TestObs);
        
    end
    for iTrial = 1:size(TestObs,2)
        Result(iTrial).fr_proj = squeeze(fr_proj(iTrial,:,:));
        Result(iTrial).fl_proj = squeeze(fl_proj(iTrial,:,:));
    end
end
if ~splitTrainTest
    save(strcat(fileHeader,'_train.mat'),'Result','EstParam','ll','InitParam','-v7.3');
end

end

%% ========== Nested functions =====================
% Function looparam removes the relevant neuron information from parameters
function LooEstParam = looparam(EstParam,mi)
LooEstParam = EstParam;
LooEstParam.A = EstParam.A(mi,:);
LooEstParam.R = EstParam.R(mi,mi);
LooEstParam.b = EstParam.b(mi);
LooEstParam.P = EstParam.P;
LooEstParam.D = EstParam.D;
LooEstParam.h_1 = EstParam.h_1;
LooEstParam.G_1 = EstParam.G_1;
end

% Function loof predicts the left out neuron's fluorescence
function [fr_proj,fl_proj] = loof(jObs,EstParam,zLoo,N_TRIAL,TestObs)
fr_proj = zeros(N_TRIAL,size(zLoo,3));
fl_proj = zeros(N_TRIAL,size(zLoo,3));

A = [EstParam.A, EstParam.b];
if isfield(TestObs,'gam')
    G =  diag(mean([TestObs(:).gam],2)); % For deconv-LDS
else
    G = nan; % for LDS direct on fluorescence
end
for nTrial = 1:N_TRIAL
    fr_proj(nTrial,:) = A(jObs,:) * squeeze(zLoo(nTrial,:,:));
    if isfield(TestObs,'gam') % deconv-LDS
        fl_proj(nTrial,1) = TestObs(1).c(jObs,1);
        for jTime = 2:size(TestObs(1).y,2)
            fl_proj(nTrial,jTime) = G(jObs,jObs)*fl_proj(nTrial,jTime-1) + ...
                A(jObs,:)*squeeze(zLoo(nTrial,:,jTime))';
        end
        fl_proj(nTrial,:) = fl_proj(nTrial,:)+TestObs(nTrial).b(jObs);
    else % LDS direct on fluorescence
        fl_proj(nTrial,:) = fr_proj(nTrial,:);
    end
end
end
