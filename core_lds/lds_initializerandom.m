function Param = lds_initializerandom(N_LATENT,N_NEURON,varargin)
%  Randomly generates the parameters D, P, A, R, b
%  of the LDS system.
%                z(t) = Dz(t-1) + u(t)
%                y(t) = Az(t) + b + w(t)
%  where w(t) ~ N(0,R), u(t) ~ N(0,P), z(1)~N(h_1,G_1)
%
%  INPUT:
%       N_LATENT - scalar indicating dimension of latents
%       N_NEURON - scalar indicating dimension of observed data
%
%  OUTPUT:
%       Param - structure containing D, P, A, b, R parameters
%
%  OPTIONAL INPUT:
%       coeffA - scalar indicating scale of A
%       coeffD - scalar indicating scale of D (limited by eigenvalues)
%       coeffR - scalar indicating scale of R (noise term higher dim)
%       coeffP - scalar indicating scale of P
%       coeffb - scalar indicating scale of b term
%% AUTHOR    : Koh Tze Hui
%% DEVELOPED : MATLAB (R2018a)
%% FILENAME  : lds_initializerandom.m
%% LAST CHECKED: 211101 (YYMMDD)
%% === Initialize function parameters ====
p = inputParser;
p.addRequired('N_LATENT',@isscalar);
p.addRequired('N_NEURON',@isscalar);
p.addParameter('coeffA',1,@isscalar);
p.addParameter('coeffR',0.5,@isscalar);
p.addParameter('coeffP',0.5,@isscalar);
p.addParameter('coeffb',5,@isscalar);
p.addParameter('coeff_h1',0,@isscalar);
p.addParameter('coeffG1',1,@isscalar)
p.addParameter('InitParam',struct,@isstruct);

p.parse(N_LATENT,N_NEURON,varargin{:});

N_LATENT = p.Results.N_LATENT;
N_NEURON = p.Results.N_NEURON;
coeffA = p.Results.coeffA;
coeffR = p.Results.coeffR;
coeffP = p.Results.coeffP;
coeffb = p.Results.coeffb;
coeffh1 = p.Results.coeff_h1;
coeffG1 = p.Results.coeffG1;
InitParam = p.Results.InitParam;


%% ==== z(t) = Dz(t-1)+ u(t) ===
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

%% === y(t) = Az(t) + b + w(t) ===
if isfield(InitParam,'A')
    A = InitParam.A;
else
    A = coeffA*rand(N_NEURON,N_LATENT)-coeffA/2;
end

if isfield(InitParam,'b')
    b = InitParam.b;
else
    b = coeffb*rand(N_NEURON,1);
end

if isfield(InitParam,'R')
    R = InitParam.R;
else
    R = coeffR*diag(rand(N_NEURON,1));
end

%% === z(1) ===
if isfield(InitParam,'G_1')
    G_1 = InitParam.G_1;
else
    G_1 = coeffG1*diag(rand(N_LATENT,1));
end
 
if isfield(InitParam,'h_1')
    h_1 = InitParam.h_1;
else
    h_1 = coeffh1*rand(N_LATENT,1);
end

Param.D = D; Param.P = P; Param.A = A; Param.R = R; Param.b = b;
Param.h_1 = h_1; Param.G_1 = G_1;
end


