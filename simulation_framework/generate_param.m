%% Function generate_param
%  Function generate_param generates the parameters G, A, b, Q, B, d, R,...
%  D, k, P, mu_1, cov_1 for the model:
%                z(t) = Dz(t-1) + k + u(t)
%                c(t) = Gc(t-1) + Az(t) + b + v(t)
%                y(t) = Bc(t) + d + w(t)
%  where v(t) ~ N(0,Q), w(t) ~ N(0,R), u(t) ~ N(0,P)
%  c(1) ~ N(mu_1,cov_1), z(1)~N(0,I)
%
%  Input:
%       P - scalar indicating dimension of latents
%       D - scalar indicating dimension of observed data      
%  
%  Parameters:
%       coeffB - scalar indicating scale of B (diag)
%       coeffA - scalar indicating scale of A  
%       coeffD - scalar indicating scale of D (diag)
%       coeffR - scalar indicating scale of R (diag)
%       coeffQ - scalar indicating scale of Q (diag)
%       coeffP - scalar indicating scale of P (diag)
%       coeffMu - scalar indicating scale of mu_1
%       coeffCov - scalar indicating scale of cov_1
%       coeff_d - scalar indicating scale of bias vector d
%       coeff_b - scalar indicating scale of bias vector b   
%       coeff_k - scalar indicating scale of bias vector k
%       FixParam - struct indicating fixed parameters
%
%  Output:
%       Param - structure containing G, A, b, Q, B, d, R, D, k, P, ...
%               cov_1, mu_1 parameters        
%% AUTHOR    : Koh Tze Hui
%% DEVELOPED : MATLAB (R2018a)
%% FILENAME  : generate_param.m

function Param = generate_param(N_latent,N_proj,varargin)
p = inputParser;
p.addRequired('N_latent',@isscalar);
p.addRequired('N_proj',@isscalar);
p.addParameter('coeffB',1,@isscalar);
p.addParameter('coeffR',1000,@isscalar);
p.addParameter('coeffA',5,@isscalar);
p.addParameter('coeffQ',1,@isscalar);
p.addParameter('coeffD',1,@isscalar);
p.addParameter('coeffP',1,@isscalar);
p.addParameter('coeff_b',0.1,@isscalar);
p.addParameter('shift_b',0,@isscalar);
p.addParameter('coeff_d',1,@isscalar);
p.addParameter('coeff_k',1,@isscalar);
p.addParameter('coeffMu',5,@isscalar);
p.addParameter('coeffCov',5,@isscalar);
p.addParameter('InitParam',struct,@isstruct);
p.parse(N_latent,N_proj,varargin{:});

N_latent = p.Results.N_latent;
N_proj = p.Results.N_proj;
coeffB = p.Results.coeffB;
coeffR = p.Results.coeffR;
coeffA = p.Results.coeffA;
coeffQ = p.Results.coeffQ;
coeffD = p.Results.coeffD;
coeffP = p.Results.coeffP;
coeffb = p.Results.coeff_b;
shiftb = p.Results.shift_b;
coeffd = p.Results.coeff_d;
coeffk = p.Results.coeff_k;
coeffMu = p.Results.coeffMu;
coeffCov = p.Results.coeffCov;
InitParam = p.Results.InitParam;

% Assign parameters if fixed

%% y(t) = Bc(t) + d + w(t)
if isfield(InitParam,'B')
    B = InitParam.B;
else
    B = coeffB*diag(rand(N_proj,1));
end

if isfield(InitParam,'d')
    d = InitParam.d;
else
    d = coeffd*rand(N_proj,1);
end

if isfield(InitParam,'R')
    R = InitParam.R;
else
    R = coeffR*diag(rand(N_proj,1));
end
%% c(t) = Gc(t-1) + Az(t) + b + v(t)
if isfield(InitParam,'G')
    G = InitParam.G;
else
    G = diag(0.005*rand(N_proj,1) + 0.995); % Range 0.5-1.0
end

if isfield(InitParam, 'A')
    A = InitParam.A;
else
    A = generate_a(N_proj,N_latent,coeffA);
end

if isfield(InitParam,'b')
    b = InitParam.b;
else
    b = coeffb*rand(N_proj,1)+shiftb;
end

if isfield(InitParam,'Q')
    Q = InitParam.Q;
else
    Q = coeffQ*diag(rand(N_proj,1));
end

%% z(t) = Dz(t-1) + k + u(t)
if isfield(InitParam,'D')
    D = InitParam.D;
else
%     D = coeffD*diag(rand(N_latent,1));
    D = diag(0.005*rand(N_latent,1) + 0.995); % Range 0.5-1.0
end

if isfield(InitParam,'P')
    P = InitParam.P;
else
    P = coeffP*diag(rand(N_latent,1));
end

if isfield(InitParam,'k')
    k = InitParam.k;
else
    k = coeffk*randn(N_latent,1);
end

%% c(1)
if isfield(InitParam,'cov_1')
    cov_1 = InitParam.cov_1;
else
    cov_1 = generate_cov(N_proj, coeffCov);
end
 
if isfield(InitParam,'mu_1')
    mu_1 = InitParam.mu_1;
else
    mu_1 = coeffMu*rand(N_proj,1);
%     mu_1 = zeros(N_proj,1);
end



Param.B = B; Param.R = R; Param.d = d;
Param.G = G; Param.A = A; Param.Q = Q;  Param.b = b;
Param.D = D; Param.P = P; Param.k = k;
Param.cov_1 = cov_1; Param.mu_1 = mu_1; 
end