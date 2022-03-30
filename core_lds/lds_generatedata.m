%% Function lds_generatedata
%  Function lds_generatedata generates the latents (z) and observations (y) of
%  the LDS system using the specified parameters according to the equations
%           z(t+1) = D*z(t) + u(t)
%           y(t)   = A*z(t) + b + w(t)
%  where u(t)~N(0,P), w(t)~N(0,R), z(1)~N(h_1,G_1)
%
%  Input:
%       GenParam - Structure containing parameters D, A, P, R, b
%  
%  Output:
%       Observation - 1xN_TRIAL structure containing y, the observed data
%       Latent      - 1xN_TRIAL structure containing z, the lataent data
%
%  Parameters:
%       Latent - Structure containing generated latents
%       T     - scalar indicating length of trial
%       N_TRIAL     - scalar indicating number of trials
%
%% AUTHOR    : Koh Tze Hui
%% DEVELOPED : MATLAB (R2018a)
%% FILENAME  : lds_generate_data.m
%% LAST CHECKED: 200420 (YYMMDD)

function [Observation, Latent] = lds_generatedata(GenParam,varargin)
p = inputParser;
p.addRequired('GenParam',@isstruct);
p.addParameter('Latent',struct,@isstruct);
p.addParameter('N_TRIAL',100,@isscalar);
p.addParameter('T',100,@isscalar);
p.parse(GenParam,varargin{:});
GenParam = p.Results.GenParam;
Latent = p.Results.Latent;

if size(Latent,2)==1
    N_TRIAL = p.Results.N_TRIAL;
    T = p.Results.T;
    latentInput = false;
else
    N_TRIAL = size(Latent,2);
    T = size(Latent(1).x,2);
    latentInput = true;
end
R = GenParam.R; P = GenParam.P; D = GenParam.D; A = GenParam.A; b = GenParam.b;
h_1 = GenParam.h_1; G_1 = GenParam.G_1;
N_LATENT = size(P,1); N_NEURON = size(R,1);

Observation = struct; Latent = struct;

for iTrial = 1:N_TRIAL
    y = zeros(N_NEURON,T);
    z = zeros(N_LATENT,T);
    u = mvnrnd(zeros(T,N_LATENT),P)';
    w = mvnrnd(zeros(T,N_NEURON),R)';
    if ~latentInput
        z(:,1) = mvnrnd(h_1,G_1)';
    else
        z = Latent(iTrial).z;
    end
    y(:,1) = A*z(:,1) + b + w(:,1);
    for jTime = 2:T
        z(:,jTime) = D*z(:,jTime-1) + u(:,jTime);
        y(:,jTime) = A*z(:,jTime) + b + w(:,jTime);
    end
    Latent(iTrial).z = z;
    Observation(iTrial).y = y;
end
end