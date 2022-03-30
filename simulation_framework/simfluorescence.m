%% Function SIMFLUORESCENCE
%  Simulates a fluorescence trace for all trials
%
%  Input:
%       Param  - Structure containing commonly used parameters
%       Spktrain - Structure of spike trains
%                   dim(1xN_TRIAL)
%
% Output:
%       Ftrace - Generated fluorescence trace
%                           dim(1xN_TRIAL)
%
% Parameters:
%       kernel - selected kernel to be used
%                   1) 'noNoise'
%                   2) 'poissonNoise'
%                   3) 'gaussianNoise' (default)
%
%% AUTHOR    : Koh Tze Hui
%% DEVELOPED : MATLAB (R2017a)
%% FILENAME  : simfluorescence.m

function [Ftrace] = simfluorescence(Param,Spktrain,gammaIdx,varargin)
p = inputParser;
p.addRequired('Param',@isstruct);
p.addRequired('Spktrain',@(x) isstruct(x));
p.addRequired('gammaIdx',@isscalar);
p.addParameter('kernel','gaussianNoise',@(x) ischar(x));
p.addParameter('noise',nan,@isnumeric);
p.parse(Param,Spktrain,gammaIdx,varargin{:});

Param = p.Results.Param;
Spktrain = p.Results.Spktrain;
gammaIdx = p.Results.gammaIdx;
kernel = p.Results.kernel; % Selected generation method
noise = p.Results.noise;

%% Generate N_trial of fluorescence
Ftrace = struct;
for iSplit = 1:Param.N_SPLIT
    if ~isnan(noise)
        [y,G,c] = convolvefluorescence(Param,Spktrain(iSplit).y,gammaIdx,'kernel',kernel,'noise',noise);
    else
        [y,G,c] = convolvefluorescence(Param,Spktrain(iSplit).y,gammaIdx,'kernel',kernel);
    end
    Ftrace(iSplit).y = y;
    Ftrace(iSplit).c = c;
end
Ftrace(1).G = G;
disp('All fluorescence traces have been simulated');
end