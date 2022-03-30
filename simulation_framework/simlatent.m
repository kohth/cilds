%% Function SIMLATENT
%  Simulates x,the N number of latents using varying processes E.g. Gaussian
%  process
%
% Input:
%       Param - Structure containing commonly used parameters
%       tau - timescales to be used in generating latents
%                   dim(1xN_LATENT)
%
% Output:
%       z - generated gaussian process in the form of a 1xn array
%
% Parameters:
%       kernel - String containing selected generation method
%                       1) 'gaussianprocess' (default)

%% AUTHOR    : Koh Tze Hui
%% DEVELOPED : MATLAB (R2017a)
%% FILENAME  : simlatent.m

function x = simlatent(Param,tau,varargin)
p = inputParser;
p.addRequired('Param',@isstruct);
p.addRequired('tau',@(x) isnumeric(x));
p.parse(Param,tau,varargin{:});

Param = p.Results.Param;
tau = p.Results.tau;

kernel = Param.LATENTTYPE;
x = zeros(Param.N_LATENT,Param.N_TRIAL*Param.T_GEN/Param.N_SPLIT);


switch kernel
    case 'gaussianprocess' % Gaussian process trajectory
        for i = 1:Param.N_LATENT
            x(i,:) = gaussprocess(Param,tau(i));
            disp(i);
        end
    case 'autoregressive'
        x = autoregressiveprocess(Param,tau);
    otherwise
        error('User Defined Error: Unknown latent generation kernel specified');
end


end