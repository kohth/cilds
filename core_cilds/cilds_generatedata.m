%% Function cilds_generatedata
%  Function cilds_generatedata generates the latents (c) (z) and observations 
%  (y) of the LDS system using the specified parameters according to the 
%  equations
%                z(t) = Dz(t-1) + u(t)
%                c(t) = Gc(t-1) + Az(t) + b + v(t)
%                y(t) = Bc(t) + w(t)
%  where v(t) ~ N(0,Q), w(t) ~ N(0,R), u(t) ~ N(0,P)
%  c(1) ~ N(mu_1,cov_1), z(1)~N(h_2,G_2)
%
%  Input:
%       GenParam - Structure containing parameters G, A, b, Q, B, R,... 
%               D, P, mu_1,cov_1,h_2,G_2
%
% Output:
%       Observation - 1xN_TRIAL structure containing y, the observed data
%       Latent      - 1xN_TRIAL structure containing z and c, the latent data
%
% OPTIONAL INPUT:
%       Latent   - 1xN_TRIAL Structure containing latent variables z (default: empty structure)
%       N_TRIAL  - scalar indicating number of trials (default: 100)
%       T        - scalar indicating trial length (default: 100)
%
%% AUTHOR    : Koh Tze Hui
%% DEVELOPED : MATLAB (R2018a)
%% FILENAME  : lds_generatedata.m
%% LAST CHECKED: 211101 (YYMMDD)

function [Observation, Latent] = cilds_generatedata(GenParam,varargin)
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
    T = size(Latent(1).z,2);
    latentInput = true;
end
B = GenParam.B; R = GenParam.R; 
G = GenParam.G; A = GenParam.A; b = GenParam.b; Q = GenParam.Q; 
D = GenParam.D; P = GenParam.P;
mu_1 = GenParam.mu_1; cov_1 = GenParam.cov_1;
h_2 = GenParam.h_2; G_2 = GenParam.G_2;
N_NEURON = size(B,1); N_LATENT = size(D,1);

Observation = struct; 

for iTrial = 1:N_TRIAL
    
    y = zeros(N_NEURON,T);
    c = zeros(N_NEURON,T);
    if ~latentInput
    z = zeros(N_LATENT,T);
    end
    v = mvnrnd(zeros(1,N_NEURON),Q,T)';
    w = mvnrnd(zeros(1,N_NEURON),R,T)';
    u = mvnrnd(zeros(1,N_LATENT),P,T)';
    c(:,1) = mvnrnd(mu_1,cov_1)';
    if ~latentInput
    z(:,2) = mvnrnd(h_2,G_2)'; % c_1 doesn't depend on z_1 so no need for z_1
    else
        z = Latent(iTrial).z;
    end
    y(:,1) = B*c(:,1)+ w(:,1);
    for jTime = 2:T
        if ~latentInput && jTime > 2
        z(:,jTime) = D*z(:,jTime-1) + u(:,jTime);
        end
        c(:,jTime) = G*c(:,jTime-1) + A*z(:,jTime) + b + v(:,jTime);
        y(:,jTime) = B*c(:,jTime) + w(:,jTime);
    end
    Latent(iTrial).c = c;

    Latent(iTrial).z = z; 
    Observation(iTrial).y = y;
end