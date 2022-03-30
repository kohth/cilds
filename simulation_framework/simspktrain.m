%% Function IHPOISSPROCESS
%  Generates spike times according to an inhomogeneous poisson process using
%  method from Neural Signal Processing notes (Byron Yu)
%
% Input:
%       Param  - Structure containing commonly used parameters
%       Model  - Structure containing field with of time varying
%                firing rates
%                   dim(1xN_TRIAL)
%
%
% Output:
%       Spktrain - spike train generated from an inhomogeneous poisson
%                process
%                   dim(1xN_TRIAL)
%
% Parameters:
%       units - string indicating whether lambda is in seconds or
%               milliseconds
%
%% AUTHOR    : Koh Tze Hui
%% DEVELOPED : MATLAB (R2018a)
%% FILENAME  : simspktrain.m

function Spktrain = simspktrain(Param,Model,varargin)
p = inputParser;
p.addRequired('Param',@isstruct);
p.addRequired('Model',@isstruct);
p.addParameter('units','s',@ischar);
p.parse(Param,Model,varargin{:});

Param = p.Results.Param;
Model = p.Results.Model;
units = p.Results.units;

%% Generate N_trial of spike trains

Spktrain = struct;

for iSplit = 1:Param.N_SPLIT
    Spktrain(iSplit).trialId = iSplit;
    if strcmp(units,'s')
        lambda = num2cell(Model(iSplit).y,2);
    elseif strcmp(units,'ms')
        lambda = num2cell(1000.*(Model(iSplit).y),2);
    else
        error('Unknown units');
    end
    yCurr = cellfun(@(x) ihpoissprocess(Param,x),lambda,'un',0);
    Spktrain(iSplit).y = cell2mat(yCurr);
    disp(strcat('Spike Train trial: ',num2str(iSplit)));
end

disp('All spike trains have been simulated');

end
