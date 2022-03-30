%% Function SIMOFFSET
%  Simulates d, the offset matrix
%
% Input:
%       Param - Structure containing commonly used parameters
%
% Output:
%       d - generated offset matrix 
%               dim(N_PROJ x T_GEN)
%
% Parameters:
%       kernel - String containing selected generation method
%                       1) 'gaussian' (default)
%       isTest - Logical for plotting out the generated vector
%                       1) false (default)
%% AUTHOR    : Koh Tze Hui
%% DEVELOPED : MATLAB (R2017a)
%% FILENAME  : simloading.m

function d = simoffset(Param,varargin)
%% Parameters
p = inputParser;
p.addRequired('Param',@isstruct);
p.addParameter('kernel','uniform',@(x) isstring(x)|| ischar(x));
p.addParameter('isTest',false,@islogical);
p.parse(Param,varargin{:});

Param = p.Results.Param;
kernel = p.Results.kernel; % Selected generation method
isTest = p.Results.isTest;

switch kernel
    case 'uniform'
        d = repmat(Param.DSHIFT,Param.N_PROJ,1) + Param.DVAR*rand(Param.N_PROJ,1);
    case 'gaussian'
        d = repmat(Param.DSHIFT,Param.N_PROJ,1) + ...
            Param.DVAR*randn(Param.N_PROJ,1);
    otherwise
        error('User Defined Error: Unknown offset kernel selected');
end

% For debugging
if isTest
    figure();
    plot(d);
    title('Offset vector'); xlabel('Neuron No.'); ylabel('Offset');
    figName = strcat('./',Param.PLOTFOLDER,'/offset_vector.png');
    print(figName,'-dpng');
end
end