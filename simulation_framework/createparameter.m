%% Function CREATEPARAMETER 
%  creates the generally used parameters in a structure
%
%  Input: 
%       -
%  Output: Structure
%       Param.N_TRIAL                 - Number of trials
%
%       Param.N_PROJ                  - Number of projected neurons
%
%       Param.N_LATENT                - Number of latents
%
%       Param.T_GEN                   - Length of generated trial (msec)
%
%       Param.T_MAX                   - Length to truncate trial to (msec)
%
%       Param.T_STEP                  - Scalar time step (msec)
%
%       Param.T_VEC                   - Time vector (msec) 
%                                           dim(1xT_GEN)
%
%       Param.BIN                     - Scalar bin width (msec)
%
%       Param.SEPARATETRIAL           - Boolean indicating whether short,
%                                       independent trials should be used
%                                       (fluorescence starts from 0 each
%                                       trial), or a long gaussian trace 
%                                       segmented into trials
%
%       Param.BIN_VEC                 - Vector containing time points of
%                                       bins
%                                           dim(1xT_MAX/BIN)
%
%       Param.TAU                     - Array of taus to test on
%                                           dim(1xnTau)
%
%       Param.BSHIFT                  - Mean of the offset vector
%
%       Param.BVAR                    - Variance of the offset vector
%
%       Param.ASHIFT                  - Mean of the loading matrix
%
%       Param.AVAR                    - Variance of the loading matrix
%
%       Param.GAMMA                   - Fluorescence AR gamma array
%
%       Param.GAMMAKNOWN              - Tell deconvolution if gamma is
%                                       known
%
%       Param.NOISE                   - Fluorescence gaussian noise level
%
%       Param.DECONVOLUTIONTYPE       - Selected type of deconvolution
%
%       Param.SAMPLEDDATAFOLDER       - Name of folder to store sampled
%                                       data
%
%       Param.DATAFOLDER              - Name of folder to store data
%
%       Param.PLOTFOLDER              - Name of folder to store plots
%
%       Param.RESULTFOLDER            - Name of folder to store results
%
%       Param.INITTYPE                - Initialization type for
%                                       dimensionality reduction
%
%       Param.LATENTTYPE              - Type of latents generated
%
%       Param.GENTYPE                 - Type of observations generated
%
%% AUTHOR    : Koh Tze Hui
%% DEVELOPED : MATLAB (R2018a)
%% FILENAME  : createparameter.m
    
function RunParam = createparameter()
%% Store information into a structure for output
RunParam.N_TRIAL =100; %200 in paper
RunParam.N_REMOVE = 0;
RunParam.N_PROJ = 94;
RunParam.N_LATENT = 10;
RunParam.T_GEN = 30000; %60000 in paper
RunParam.T_MAX = 29000; %59000 in paper
RunParam.T_STEP = 1;

RunParam.PIECESIZE = 5000;
RunParam.HISTORY = 5000;
RunParam.BIN = 25;
RunParam.N_SPLIT = 2;
T = RunParam.T_GEN*RunParam.N_TRIAL/RunParam.N_SPLIT;
RunParam.T_VEC = 0:RunParam.T_STEP:T-RunParam.T_STEP;
RunParam.BIN_VEC = floor((RunParam.BIN+1)/2):RunParam.BIN:T-floor((RunParam.BIN-1)/2);
RunParam.TAU = [50,100,200,1000,2000,5000];
RunParam.GAMMA = [0.9985,0.9993,0.9996];

RunParam.NOISE = 1.5;
%% Deconvolution
RunParam.GAMMAKNOWN = false;
RunParam.DECONVOLUTIONTYPE = 'constrained';
RunParam.SMIN = 0;
RunParam.OPTIMIZEPARS = true;
RunParam.OPTIMIZEB = true;
RunParam.OPTIMIZESMIN = true;
RunParam.INITIALIZEGAMMA = true;
%% File names
RunParam.SAMPLEDDATAFOLDER = 'sim_datasample';
RunParam.FULLDATAFOLDER = 'sim_datafull';
RunParam.PLOTFOLDER = 'sim_plots';
RunParam.RESULTFOLDER = 'sim_results';

RunParam.INITTYPE = 'ldsInit';
RunParam.LATENTTYPE = 'gaussianprocess';

if strcmp(RunParam.LATENTTYPE, 'autoregressive')
    RunParam.BSHIFT = 0.00;
    RunParam.BVAR = 0.01;
    RunParam.AVAR = 0.005;
    RunParam.DCONST = [0.924743408289499;0.996881832383117;0.999218873791583;0.999801886843323;0.999992688731269;0.999996520297324];
elseif strcmp(RunParam.LATENTTYPE, 'gaussianprocess')
    RunParam.BSHIFT = 10;
    RunParam.BVAR = 15;
    RunParam.AVAR = 10;
    
end
RunParam.REALDATA = "M1";
RunParam.PARTIALSAVE = true;

if RunParam.N_SPLIT<RunParam.N_TRIAL
RunParam.SPLITIND = randsample(RunParam.N_SPLIT,floor(RunParam.N_SPLIT/2));
RunParam.TRAININD = reshape(1:RunParam.N_TRIAL,[RunParam.N_TRIAL/RunParam.N_SPLIT,RunParam.N_SPLIT]);
RunParam.TRAININD = RunParam.TRAININD(:,RunParam.SPLITIND);
RunParam.TRAININD = [RunParam.TRAININD(:)];
else
RunParam.TRAININD =  randsample(RunParam.N_TRIAL,floor(RunParam.N_TRIAL/2));
RunParam.SPLITIND = RunParam.TRAININD;
end
end

