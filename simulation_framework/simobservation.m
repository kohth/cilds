%% Function SIMOBSERVATION
%  Simulates y, the firing rates, from the latents, loading and offset
%  matrices
%
% Input:
%       Param - Structure containing commonly used parameters
%       x     - Vector of latents
%                   dim(N_LATENTxT_GEN)
%       A     - Loading matrix
%                   dim(N_PROJxN_LATENT)
%       b     - Offset matrix
%                   dim(N_PROJx1)
%
% Output:
%       y - latents projected into higher dimensions to form firing rates
%               dim(N_PROJxT_GEN)
%
%
%% AUTHOR    : Koh Tze Hui
%% DEVELOPED : MATLAB (R2017a)
%% FILENAME  : simfiringrate.m

function y = simobservation(Param,A,z,b,varargin)
p = inputParser;
p.addRequired('Param',@isstruct);
p.addRequired('A',@(x) isnumeric(x));
p.addRequired('z',@(x) isnumeric(x));
p.addRequired('b',@(x) isnumeric(x));

p.parse(Param,A,z,b,varargin{:});

Param = p.Results.Param;
A = p.Results.A;
z = p.Results.z;
b = p.Results.b;

b = repmat(b,1,Param.T_GEN*Param.N_TRIAL/Param.N_SPLIT);


%% Project latents into higher dimensions

y = A*z + b;


% Change negative values to positive
y = log(1+exp(y));


end