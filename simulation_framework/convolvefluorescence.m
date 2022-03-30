%% Function CONVOLVEFLUORESCENCE
%  Simulates a fluorescence trace given a binary spike train by convolving
%  with a kernel
%
%  Input:
%       Param  - Structure containing commonly used parameters
%       spkTrain - Spike trains where N_PROJ is the number of
%                neurons and T_GEN is the length of the train
%                   dim(N_PROJxT_GEN)
%
% Output:
%       fluorescence - Generated fluorescence trace
%                           dim(N_PROJXT_GEN)
%
% Parameters:
%       kernel - selected kernel to be used
%                   1) 'noNoise'
%                   2) 'poissonNoise'
%                   3) 'gaussianNoise' (default)
%
%% AUTHOR    : Koh Tze Hui
%% DEVELOPED : MATLAB (R2017a)
%% FILENAME  : convolvefluorescence.m


function [fluorescence,G,truth] = convolvefluorescence(Param,spkTrain,gammaIdx,varargin)
p = inputParser;
p.addRequired('Param',@isstruct);
p.addRequired('spkTrain',@(x) isnumeric(x));
p.addRequired('gammaIdx',@isscalar);
p.addParameter('kernel','gaussianNoise',@(x) ischar(x));
p.addParameter('noise',nan,@isnumeric);
p.addParameter('jitter',false,@islogical);
p.parse(Param,spkTrain,gammaIdx,varargin{:});

Param = p.Results.Param;
spkTrain = p.Results.spkTrain;
gammaIdx = p.Results.gammaIdx;
noise = p.Results.noise;
jitter = p.Results.jitter;
kernel = p.Results.kernel; % Selected generation method

b = 0;
p = length(Param.GAMMA(gammaIdx));
gamma = [flipud(reshape(Param.GAMMA(gammaIdx), [], 1)); 1];

T = Param.T_GEN*Param.N_TRIAL/Param.N_SPLIT;


% for t=(p+1):T
%     truth(:, t) = truth(:, (t-p):t) * gamma;
% end


truth = double(spkTrain);
if jitter
    G = diag(Param.GAMMA(gammaIdx)*ones(Param.N_PROJ,1)+(0.0005*rand(Param.N_PROJ,1)-0.0005));
else
    G = diag(Param.GAMMA(gammaIdx)*ones(Param.N_PROJ,1));
end
for t=(p+1):T
    truth(:, t) = G*truth(:,t-1)+truth(:,t);
end


if ~isnan(noise)
    coeff = noise;
else
    coeff = Param.NOISE;
end

switch kernel
    case 'noNoise' % Clean exponential decay
        
    case 'poissonNoise' % poisson noise added
        % Define parameters according to theis,2016
        a = 100;
        fluorescence = poissrnd(b + a*truth);
        
    case 'gaussianNoise' % gaussian noise added
        % Define parameters according to vogelstein,2010
        fluorescence = b + truth + coeff * randn(Param.N_PROJ, T);
        
end

end