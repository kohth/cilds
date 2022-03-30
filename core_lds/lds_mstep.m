function EstParam =lds_mstep(Xf,Vf,Vj,PrevParam,Observation,FixParam)
% Perform maximization step of EM algorithm
%
% INPUT:
%
%      Xf           - 3D matrix containing E[z_t|y_1,...,y_T] for all T,
%                           - Dimensions: N_TRIAL x (N_LATENT) x T
%
%      Vf           - 3D matrix containing cov[z_t|y_1,...,y_T]
%                           - Dimensions: N_LATENT x N_LATENT x T
%
%      Vj           - 3D matrix containing cov[z_t,z_(t-1)}y_1,....,y_T]
%                           - Dimensions: N_LATENT x N_LATENT x T
%
%      PrevParam    - Structure containing model parameter from previous
%                     maximization iteration, since we are doing partial
%                     maximization each time
%                           - Dimensions: 1 x 1
%                           - Fields: 
%                               A (N_NEURONxN_LATENT),
%                               D (N_LATENTxN_LATENT),
%                               R (N_NEURONxN_NEURON),
%                               P (N_NEURONxN_LATENT),
%                               b (N_NEURONx1),
%                               h_1 (N_LATENTx1),
%                               G_1 (N_LATENTxN_LATENT)
%
%      Observation  - Structure containing recorded data (fluorescence
%                     traces)
%                           -Dimensions: 1 x N_TRIAL
%                           -Fields:
%                               y (N_NEURON x T) -- neural data
%
%      FixParam     - Structure indicating which parameters to hold
%                     constant in M step (i.e. don't update the values
%                     beyond initialization)
%                           - Default: empty
%                           - Possible fields:
%                               A,D,R,P,b,h_1,G_1
% 
%  OUTPUT:
%      EstParam    - Structure containing model parameter after current
%                    maximization step
%                           - Dimensions: 1 x 1
%                           - Fields:
%                               A,D,R,P,b,h_1,G_1
%
%% AUTHOR    : Koh Tze Hui
%% DEVELOPED : MATLAB (R2018a)
%% FILENAME  : lds_mstep.m
%% LAST CHECKED: 211101 (YYMMDD)
%% REFERENCES:
%     Ghahramani Z Parameter Estimation for Linear Dynamical Systems

NLatent = size(Xf,2); NProj= size(Observation(1).y,1);
NTrial = size(Observation,2);
T = size(Observation(1).y,2);

%% === Initialization for summations ===
sumyx = zeros(NProj,NLatent);
sumPt = zeros(NLatent,NLatent);
sumPtt_1 = zeros(NLatent,NLatent);
sumyy = zeros(NProj,NProj);
y = zeros(NTrial,NProj,size(Observation(1).y,2));
for iTrial = 1:NTrial
    y(iTrial,:,:) = Observation(iTrial).y;
end

%% === Summations ===
for t = T(1):-1:1
    sumPt =sumPt + Vf(:,:,t) + Xf(:,:,t)'*Xf(:,:,t)/NTrial;
    sumyx = sumyx + y(:,:,t)'*Xf(:,:,t)/NTrial;
    sumyy = sumyy + y(:,:,t)'*y(:,:,t)/NTrial;
    if t > 1
        sumPtt_1 = sumPtt_1 + Vj(:,:,t) + Xf(:,:,t)'*Xf(:,:,t-1)/NTrial;
    end
end
sumPt_1 = sumPt - Vf(:,:,T) - Xf(:,:,T)'*Xf(:,:,T)/NTrial;
sumPt_2 = sumPt - Vf(:,:,1) - Xf(:,:,1)'*Xf(:,:,1)/NTrial;

sumPt_1 = 0.5*(sumPt_1+sumPt_1');
sumPt_2 = 0.5*(sumPt_2+sumPt_2');
sumPt = 0.5*(sumPt+sumPt');
sumyy = 0.5*(sumyy+sumyy');

%% A, R, D, P
if ~isfield(FixParam,'D')
    D = diag(diag(sumPtt_1))/(diag(diag(sumPt_1)));
else
    D = PrevParam.D;
end
if ~isfield(FixParam,'A')
    A = sumyx/sumPt;
else
    A = PrevParam.A;
end

if ~isfield(FixParam,'R')
    R = (1/T)*(diag(diag(sumyy)) - diag(diag(A*sumyx')));
else
    R = PrevParam.R;
end
if ~isfield(FixParam,'P')
    P = (1/(T-1))*(diag(diag(sumPt_2)) - diag(diag(D* sumPtt_1')));
else
    P = PrevParam.P;
end

%% mu_1, cov_1
if ~isfield(FixParam,'h_1')
    h_1 = (sum(Xf(:,:,1),1))';
    h_1 = h_1/NTrial;
else
    h_1 = FixParam.h_1;
end

if ~isfield(FixParam,'G_1')
    G_1 = diag(diag(Vf(:,:,1))) + diag(diag((Xf(:,:,1)-repmat(h_1,1,NTrial)')'*(Xf(:,:,1)-repmat(h_1,1,NTrial)')))/NTrial;
else
    G_1 = FixParam.G_1;
end
G_1 = 0.5*(G_1+G_1');
EstParam.G_1 = G_1;
EstParam.h_1 = h_1;


NLatent = NLatent-1;
EstParam.h_1 = h_1(1:NLatent);
EstParam.G_1 = G_1(1:NLatent,1:NLatent);
EstParam.A = A(:,1:NLatent);
EstParam.b = A(:,end);
EstParam.D = D(1:NLatent,1:NLatent);
EstParam.P = P(1:NLatent,1:NLatent);
EstParam.R = R;
end