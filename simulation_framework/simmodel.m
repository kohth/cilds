%% Function simmodel
%  Using different latents for each trial, generate firing rates for each
%  neuron of each trial
%
% Input:
%       Param - Structure containing simulation parameters
%       tau - Vector containing laent timescales
%       GenParam - Structure containing model parameters
%
% Output:
%       Model - Structure containing latents, firing rates, loading matrix
%             and offset matrix
%                 dim(1x1000)
%             Fields:
%                 1) A - Loading matrix
%                           dim(N_PROJxN_LATENT)
%                 2) b - Offset matrix
%                           dim(N_PROJx1)
%                 3) z - Latents
%                           dim(N_LATENTxT_GEN)
%                 4) y - Firing rates
%                           dim(N_PROJxT_GEN)
%                 5) trialId - Scalar identification number for the trial
%
%% AUTHOR    : Koh Tze Hui
%% DEVELOPED : MATLAB (R2018a)
%% FILENAME  : simmodel.m

function Model = simmodel(Param,tau,GenParam,varargin)
p = inputParser;
p.addRequired('Param',@isstruct);
p.addRequired('tau',@(x) isnumeric(x));
p.addRequired('GenParam',@isstruct);
p.parse(Param,tau,GenParam,varargin{:});

Param = p.Results.Param;
tau = p.Results.tau;
GenParam = p.Results.GenParam;

% Begin simulation of model
Model = struct;
for iSplit = 1:Param.N_SPLIT
    z = simlatent(Param,tau);
    Model(iSplit).y = simobservation(Param,GenParam.A,z,GenParam.b);
    Model(iSplit).z = z;
    Model(iSplit).trialId = iSplit;
    disp(strcat('Model trial: ',num2str(iSplit)));
end

Model(1).A = GenParam.A; Model(1).b = GenParam.b;
end

