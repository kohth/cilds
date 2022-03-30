function Param = cilds_initializerandom(N_LATENT,N_NEURON,varargin)
%  Function cilds_initializerandom randomly generates the parameters G, A, b, Q, B, R,...
%  D, P, mu_1, cov_1, h_2, G_2 for the model:
%                z(t) = Dz(t-1) + u(t)
%                c(t) = Gc(t-1) + Az(t) + b + v(t)
%                y(t) = Bc(t) + w(t)
%  where v(t) ~ N(0,Q), w(t) ~ N(0,R), u(t) ~ N(0,P)
%  c(1) ~ N(mu_1,cov_1), z(1)~N(h_2,G_2)
%
%  INPUT:
%       N_LATENT - scalar indicating dimension of latents
%       N_NEURON - scalar indicating dimension of observed data      
%
%  OUTPUT:
%       Param - structure containing G, A, b, Q, B, R, D, P, ...
%               cov_1, mu_1, G_2, h_2
%
%  OPTIONAL INPUT:
%       coeffB - scalar indicating scale of B (diag)
%       coeffA - scalar indicating scale of A  
%       coeffD - scalar indicating scale of D (diag)
%       coeffR - scalar indicating scale of R (diag)
%       coeffQ - scalar indicating scale of Q (diag)
%       coeffP - scalar indicating scale of P (diag)
%       coeffMu - scalar indicating scale of mu_1
%       coeffCov - scalar indicating scale of cov_1
%       coeff_b - scalar indicating scale of bias vector b   
%       InitParam - struct indicating fixed initialal parameters
%       
%% AUTHOR    : Koh Tze Hui
%% DEVELOPED : MATLAB (R2018a)
%% FILENAME  : cilds_initializerandom.m
%% LAST CHECKED: 220322 (YYMMDD)
%% === Initialize function parameters ====

p = inputParser;
p.addRequired('N_LATENT',@isscalar);
p.addRequired('N_NEURON',@isscalar);
p.addParameter('coeffB',1,@isscalar);
p.addParameter('coeffR',0.01,@isscalar);
p.addParameter('coeffA',10,@isscalar);
p.addParameter('coeffQ',0.01,@isscalar);
p.addParameter('coeffP',0.01,@isscalar);
p.addParameter('coeff_b',1,@isscalar);
p.addParameter('coeff_mu',1,@isscalar);
p.addParameter('coeffCov',1,@isscalar);
p.addParameter('coeff_h2',0,@isscalar);
p.addParameter('coeffG2',1,@isscalar);
p.addParameter('InitParam',struct,@isstruct);
p.parse(N_LATENT,N_NEURON,varargin{:});

N_LATENT = p.Results.N_LATENT;
N_NEURON = p.Results.N_NEURON;
coeffB = p.Results.coeffB;
coeffR = p.Results.coeffR;
coeffA = p.Results.coeffA;
coeffQ = p.Results.coeffQ;
coeffP = p.Results.coeffP;
coeffb = p.Results.coeff_b;
coeffMu = p.Results.coeff_mu;
coeffCov = p.Results.coeffCov;
coeffh2 = p.Results.coeff_h2;
coeffG2 = p.Results.coeffG2;
InitParam = p.Results.InitParam;

%% y(t) = Bc(t) + w(t)
if isfield(InitParam,'B')
    B = InitParam.B;
else
    B = coeffB*diag(rand(N_NEURON,1));
end

if isfield(InitParam,'R')
    R = InitParam.R;
else
    R = coeffR*diag(rand(N_NEURON,1));
end
%% c(t) = Gc(t-1) + Az(t) + b + v(t)
if isfield(InitParam,'G')
    G = InitParam.G;
else
    G = diag(0.05*rand(N_NEURON,1) + 0.95); % Range 0.95-1.0
end

if isfield(InitParam, 'A')
    A = InitParam.A;
else
    A = coeffA*randn(N_NEURON,N_LATENT);
end

if isfield(InitParam,'b')
    b = InitParam.b;
else
    b = coeffb*rand(N_NEURON,1);
end

if isfield(InitParam,'Q')
    Q = InitParam.Q;
else
    Q = coeffQ*diag(rand(N_NEURON,1));
end

%% z(t) = Dz(t-1) + u(t)
if isfield(InitParam,'D')
    D = InitParam.D;
else
    D = diag(0.1*rand(N_LATENT,1) + 0.90);
end

if isfield(InitParam,'P')
    P = InitParam.P;
else
    P = coeffP*diag(rand(N_LATENT,1));
end

%% c(1)
if isfield(InitParam,'cov_1')
    cov_1 = InitParam.cov_1;
else
%     cov_1 = generate_cov(N_NEURON, coeffCov);
    cov_1 = coeffCov*diag(rand(N_NEURON,1));
end
 
if isfield(InitParam,'mu_1')
    mu_1 = InitParam.mu_1;
else
    mu_1 = coeffMu*rand(N_NEURON,1);
end

%% z(2)
if isfield(InitParam,'G_2')
    G_2 = InitParam.G_2;
else
%     G_2 = generate_cov(N_LATENT, coeffG2);
    G_2 = coeffG2*diag(rand(N_LATENT,1));
end
 
if isfield(InitParam,'h_2')
    h_2 = InitParam.h_2;
else
    h_2 = coeffh2*rand(N_LATENT,1);
end

Param.B = B; Param.R = R;
Param.G=G; Param.A = A; Param.Q = Q;  Param.b = b;
Param.D = D; Param.P = P;
Param.cov_1 = cov_1; Param.mu_1 = mu_1; 
Param.G_2 = G_2; Param.h_2 = h_2;
end

%% Function generate_cov
%  Function generate_cov generates a symmetric positive definite matrix
%
%  Input:
%       dim - scalar indicating dimension of covariance matrix
%       noise - scalar indicating scale of covariance matrix
%
%  Output:
%       V - dimxdim symmetric positive definite matrix

function V = generate_cov(dim,noise)
%     V = noise*randn(dim);
%     V = (V*V')/2;
% %     [V,D,W] = eig(A);
%     [~,p] = chol(V);
%     if p > 0
%         warning("Matrix is not positive definite");
%     end
    
    Q = orth(randn(dim));

    lambda = rand(dim,1);
    lambda(lambda < .1) = .1;
    lambda = lambda*noise;
    V = Q*diag(lambda)*Q';
    [~,p] = chol(V);
    if p > 0
        warning('Matrix is not positive definite');
    end
end