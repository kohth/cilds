%% Function generate_data
%  Function generate_data generates the latents (c) (z) and observations 
%  (y) of the LDS system using the specified parameters according to the 
%  equations
%                z(t) = Dz(t-1) + k + u(t)
%                c(t) = Gc(t-1) + Az(t) + b + v(t)
%                y(t) = Bc(t) + d + w(t)
%  where v(t) ~ N(0,Q), w(t) ~ N(0,R), u(t) ~ N(0,P)
%  c(1) ~ N(mu_1,cov_1), z(1)~N(0,I)
%
%  Input:
%       Param - Structure containing parameters G, A, b, Q, B, d, R,... 
%               D, k, P, mu_1,cov_1
%       T     - scalar indicating length of trial
%       N     - scalar indicating number of trials
%       isTest - boolean indicating whether to plot out latents and
%                observed data
%
% Output:
%       Observation - 1xN structure containing y, the observed data
%       Latent      - 1xN structure containing x, the lataent data
%% AUTHOR    : Koh Tze Hui
%% DEVELOPED : MATLAB (R2018a)
%% FILENAME  : generate_data.m


function [Observation, Latent] = generate_data(GenParam,varargin)
p = inputParser;
p.addRequired('GenParam',@isstruct);
p.addParameter('D',[1,0;0,1],@ismatrix);
p.addParameter('Latent',struct,@isstruct);
p.addParameter('NTrial',100,@isscalar);
p.addParameter('T',100,@isscalar);
p.parse(GenParam,varargin{:});

GenParam = p.Results.GenParam;
D = p.Results.D;
Latent = p.Results.Latent;

if ~isequal(D,[1,0;0,1])
    GenParam.D = D;
end
if size(Latent,2)==1
    NTrial = p.Results.NTrial;
    T = p.Results.T;
    latentInput = false;
else
    NTrial = size(Latent,2);
    T = size(Latent(1).z,2);
    latentInput = true;
end
B = GenParam.B; d = GenParam.d; R = GenParam.R; 
G = GenParam.G; A = GenParam.A; b = GenParam.b; Q = GenParam.Q; 
D = GenParam.D; k = GenParam.k; P = GenParam.P;
NProj = size(B,1); NLatent = size(D,1);

Observation = struct; 

for iTrial = 1:NTrial
    mu_1 = GenParam.mu_1; cov_1 = GenParam.cov_1;
    y = zeros(NProj,T);
    c = zeros(NProj,T);
    if ~latentInput
    z = zeros(NLatent,T);
    end
    v = mvnrnd(zeros(T,NProj),Q)';
    w = mvnrnd(zeros(T,NProj),R)';
    u = mvnrnd(zeros(T,NLatent),P)';
    c(:,1) = mvnrnd(mu_1,cov_1)';
    if ~latentInput
    z(:,1) = randn(NLatent,1);
    else
        z = Latent(iTrial).z;
    end
    y(:,1) = B*c(:,1) + d + w(:,1);
    for jTime = 2:T
        if ~latentInput
        z(:,jTime) = D*z(:,jTime-1) + k + u(:,jTime);
        end
        c(:,jTime) = G*c(:,jTime-1) + A*z(:,jTime) + b + v(:,jTime);
        y(:,jTime) = B*c(:,jTime) + d + w(:,jTime);
    end
    Latent(iTrial).c = c;

    Latent(iTrial).z = z;
    Observation(iTrial).y = y;
end