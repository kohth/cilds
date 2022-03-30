%% Function SIMLOADING
%  Simulates C, the loading matrix
%
% Input:
%       Param - Structure containing commonly used parameters
%
% Output:
%       c - generated loading matrix 
%               dim(N_PROJ x N_LATENT)
%
% Parameters:
%       kernel - String containing selected generation method
%                       1) 'gaussian' (default)
%       isTest - Logical for plotting out the generated matrix
%                       1) false (default)
%% AUTHOR    : Koh Tze Hui
%% DEVELOPED : MATLAB (R2017a)
%% FILENAME  : simloading.m

function c = simloading(Param,varargin)
%% Parameters
p = inputParser;
p.addRequired('Param',@isstruct);
p.addParameter('kernel','gaussian',@(x) isstring(x)|| ischar(x));
p.addParameter('isTest',false,@islogical);
p.parse(Param,varargin{:});

Param = p.Results.Param;
kernel = p.Results.kernel; % Selected generation method
isTest = p.Results.isTest; 

switch kernel
    case 'gaussian'
        % C ~ N(CSHIFT, CVAR*I)
        c = repmat(Param.CSHIFT, Param.N_PROJ,Param.N_LATENT) + ...
            Param.CVAR*randn(Param.N_PROJ,Param.N_LATENT); 
    otherwise
        error('User Defined Error: Unknown loading kernel selected');
end

% For debugging
if isTest
    figure();
    plot(c);
    title('Loading Matrix'); xlabel('Neuron No.'); ylabel('Load');
    figName = strcat('./',Param.PLOTFOLDER,'/loading_matrix.png');
    saveas(gcf,figName);
end
end