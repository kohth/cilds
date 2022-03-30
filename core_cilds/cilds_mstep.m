function [EstParam,nextChangeParam] = cilds_mstep(Xf,Vf,Vj,PrevParam,Observation,FixParam,changeParam)
% Perform maximization step of EM algorithm
% Refer to section 2.1 of supplemental methods for CILDS EM
%
% INPUT:
%
%      Xf           - 3D matrix containing E[l_t|y_1,...,y_T] for all T,
%                     where l_t = [c_t; z_(t+1)]
%                           - Dimensions: N_TRIAL x (N_LATENT + N_NEURON) x T
%
%      Vf           - 3D matrix containing cov[l_t|y_1,...,y_T]
%                           - Dimensions: (N_LATENT + N_NEURON) x (N_LATENT +
%                                         N_NEURON) x T
%
%      Vj           - 3D matrix containing cov[l_t,l_(t-1)}y_1,....,y_T]
%                           - Dimensions: (N_LATENT + N_NEURON) x (N_LATENT +
%                                         N_NEURON) x T
%
%      PrevParam    - Structure containing model parameter from previous
%                     maximization iteration, since we are doing partial
%                     maximization each time
%                           - Dimensions: 1 x 1
%                           - Fields: 
%                               A (N_NEURONxN_LATENT),
%                               B (N_NEURONxN_NEURON),
%                               G (N_NEURONxN_NEURON) - Gamma,
%                               D (N_LATENTxN_LATENT),
%                               Q (N_NEURONxN_NEURON),
%                               R (N_NEURONxN_NEURON),
%                               P (N_NEURONxN_LATENT),
%                               b (N_NEURONx1),
%                               mu_1 (N_NEURONx1),
%                               cov_1 (N_NEURONxN_NEURON)
%                               h_2 (N_LATENTx1),
%                               G_2 (N_LATENTxN_LATENT)
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
%                               A,B,G,D,Q,R,P,b,mu_1,cov_1,h_2,G_2
% 
%      changeParam - String indicating which parameters to maximize in this
%                    partial maximization step
%
% OUTPUT:
%      EstParam    - Structure containing model parameter after current
%                    maximization step
%                           - Dimensions: 1 x 1
%                           - Fields:
%                               A,B,G,D,Q,R,P,b,mu_1,cov_1,h_2,G_2
%
%      nextChangeParam - String indicating parameters to maximize in next
%                        partial maximization step
%                           - Groupings: 'BGD','A','b'
%% AUTHOR    : Koh Tze Hui
%% DEVELOPED : MATLAB (R2018a)
%% FILENAME  : cilds_mstep.m
%% LAST CHECKED: 211030 (YYMMDD)
%% REFERENCES:
%     Ghahramani Z Parameter Estimation for Linear Dynamical Systems

N_LATENT = size(PrevParam.D,1); N_NEURON = size(Observation(1).y,1);
N_TRIAL = size(Observation,2);
T = size(Observation(1).y,2);


%% Initialization for summations
sumyc = zeros(N_NEURON,N_NEURON);
sumyy = zeros(N_NEURON,N_NEURON);
sumcc_1 = zeros(N_NEURON,N_NEURON);
sumzc_1 = zeros(N_LATENT,N_NEURON);
sumcz = zeros(N_NEURON,N_LATENT);
sumz1z1 = zeros(N_LATENT,N_LATENT);
sumcc = zeros(N_NEURON,N_NEURON);
sumz1z = zeros(N_LATENT,N_LATENT);
cV = Vf(1:N_NEURON,1:N_NEURON,1:end);
zV = Vf(N_NEURON+1:N_NEURON+N_LATENT,N_NEURON+1:N_NEURON+N_LATENT,:);
cc_1 = Vj(1:N_NEURON,1:N_NEURON,:);
cz1 = Vf(1:N_NEURON,N_NEURON+1:N_NEURON+N_LATENT,:);
cz = Vj(1:N_NEURON,N_NEURON+1:N_NEURON+N_LATENT,:);
z1z = Vj(N_NEURON+1:N_NEURON+N_LATENT,N_NEURON+1:N_NEURON+N_LATENT,:);
cfut = Xf(:,1:N_NEURON,:);
z1fut = Xf(:,N_NEURON+1:N_NEURON+N_LATENT,:);

y = zeros(N_TRIAL,N_NEURON,size(Observation(1).y,2));
for iTrial = 1:N_TRIAL
    y(iTrial,:,:) = Observation(iTrial).y;
end
%% Summations

for t = T(1):-1:1
    sumcc =sumcc + cV(:,:,t) + cfut(:,:,t)'*cfut(:,:,t)/N_TRIAL;
    sumyc = sumyc + y(:,:,t)'*cfut(:,:,t)/N_TRIAL;
    sumyy = sumyy + y(:,:,t)'*y(:,:,t)/N_TRIAL;
    sumz1z1 = sumz1z1+ zV(:,:,t) + z1fut(:,:,t)'*z1fut(:,:,t)/N_TRIAL;
    if t > 1
        sumcc_1 = sumcc_1 + cc_1(:,:,t) + cfut(:,:,t)'*cfut(:,:,t-1)/N_TRIAL;
        sumzc_1 = sumzc_1 + cz1(:,:,t-1)' + z1fut(:,:,t-1)'*cfut(:,:,t-1)/N_TRIAL;
        sumcz = sumcz + cz(:,:,t) + cfut(:,:,t)'*z1fut(:,:,t-1)/N_TRIAL;
        sumz1z = sumz1z  +z1z(:,:,t) + z1fut(:,:,t)'*z1fut(:,:,t-1)/N_TRIAL;
    end
end

sumc_1c_1 = sumcc - cV(:,:,T) - cfut(:,:,T)'*cfut(:,:,T)/N_TRIAL;
sumzz = sumz1z1 - zV(:,:,T) - z1fut(:,:,T)'*z1fut(:,:,T)/N_TRIAL;
sumz1z1 = sumz1z1 - zV(:,:,1) - z1fut(:,:,1)'*z1fut(:,:,1)/N_TRIAL;

sumc_1c_1 = 0.5*(sumc_1c_1+sumc_1c_1');
sumzz = 0.5*(sumzz+sumzz');
sumz1z1 =0.5*(sumz1z1+sumz1z1');
sumcc = 0.5*(sumcc+sumcc');


sumc = sum(sum(cfut,1)/N_TRIAL,3)';
sumc_1 = sum(sum(cfut(:,:,1:end-1),1)/N_TRIAL,3)';
sumz = sum(sum(z1fut(:,:,1:end-1),1)/N_TRIAL,3)';
sumz1 = sum(sum(z1fut(:,:,2:end),1)/N_TRIAL,3)';
sumy = sum(sum(y,1)/N_TRIAL,3)';
%% B, R

if isfield(FixParam,'B')
    B = FixParam.B;
else
    B = PrevParam.B;
end


if ~isfield(FixParam,'B')
    if strcmp(changeParam,'BGD')
        B = diag(diag(sumyc))/(diag(diag(sumcc)));
    end
end

EstParam.B = B;

if ~isfield(FixParam,'R')
    R = (1/T)*(diag(diag(sumyy))-2*diag(diag(B*sumyc'))...
        +diag(diag(B*sumcc*B')));
else
    R = FixParam.R;
end
EstParam.R = R;

%% G, A, b, Q

if isfield(FixParam,'G')
    G = FixParam.G;
else
    G = PrevParam.G;
end

if isfield(FixParam,'b')
    b = FixParam.b;
else
    b = PrevParam.b;
end

if isfield(FixParam,'A')
    A = FixParam.A;
else
    A = PrevParam.A;
end


if ~isfield(FixParam,'G')
    if strcmp(changeParam,'BGD')
        G = (diag(diag(sumcc_1)) - diag(diag(A*sumzc_1)) - diag(diag(b*sumc_1')))/(diag(diag(sumc_1c_1)));
    end
end

if ~isfield(FixParam,'A')
    if strcmp(changeParam,'A')
        A = (sumcz - G*sumzc_1' - b*sumz')/(sumzz);
    end
end

if ~isfield(FixParam,'b')
    if strcmp(changeParam,'b')
        b = (1/(T-1))*(sumc -(sum(cfut(:,:,1),1)/N_TRIAL)' - G*sumc_1 - A*sumz);
    end
end

EstParam.b = b;
EstParam.A = A;
EstParam.G = G;

if ~isfield(FixParam,'Q')
    Q = (1/(T-1))*(diag(diag(sumcc - cV(:,:,1) - cfut(:,:,1)'*cfut(:,:,1)/N_TRIAL))-2*diag(diag(G*sumcc_1'))...
        -2*diag(diag(A*sumcz'))-2*diag(diag(b*(sumc-(sum(cfut(:,:,1),1)/N_TRIAL)')'))+2*diag(diag(A*sumzc_1*G'))...
        +diag(diag(G*sumc_1c_1*G'))+2*diag(diag(b*sumc_1'*G'))+diag(diag(A*sumzz*A'))...
        +2*diag(diag(b*sumz'*A'))+(T-1)*diag(diag(b*b')));
else
    Q = FixParam.Q;
end
EstParam.Q = Q;

%% D, P
if isfield(FixParam,'D')
    D = FixParam.D;
else
    D = PrevParam.D;
end

if ~isfield(FixParam,'D')
    if strcmp(changeParam,'BGD')
        D = diag(diag(sumz1z))/(diag(diag(sumzz)));
    end
end

EstParam.D = D;

if ~isfield(FixParam,'P')
    P = (1/(T-1))*(diag(diag(sumz1z1))-2*diag(diag(D*sumz1z'))...
        +diag(diag(D*sumzz*D')));
else
    P = FixParam.P;
end
EstParam.P = P;

%% mu_1, cov_1
if ~isfield(FixParam,'mu_1')
    mu_1 = (sum(cfut(:,:,1),1))';
    mu_1 = mu_1/N_TRIAL;
else
    mu_1 = FixParam.mu_1;
end

if ~isfield(FixParam,'cov_1')
    cov_1 = diag(diag(cV(:,:,1))) + diag(diag((cfut(:,:,1)-repmat(mu_1,1,N_TRIAL)')'*(cfut(:,:,1)-repmat(mu_1,1,N_TRIAL)')))/N_TRIAL;
else
    cov_1 = FixParam.cov_1;
end
cov_1 = 0.5*(cov_1+cov_1');
EstParam.cov_1 = cov_1;
EstParam.mu_1 = mu_1;


if ~isfield(FixParam,'h_2')
    h_2 = (sum(z1fut(:,:,1),1))';
    h_2 = h_2/N_TRIAL;
else
    h_2 = FixParam.h_2;
end

if ~isfield(FixParam,'G_2')
    G_2 = diag(diag(zV(:,:,1))) + diag(diag((z1fut(:,:,1)-repmat(h_2,1,N_TRIAL)')'*(z1fut(:,:,1)-repmat(h_2,1,N_TRIAL)')))/N_TRIAL;
else
    G_2 = FixParam.G_2;
end
G_2 = 0.5*(G_2+G_2');
EstParam.G_2 = G_2;
EstParam.h_2 = h_2;


if strcmp(changeParam,'BGD')
    nextChangeParam = 'A';
elseif strcmp(changeParam,'A')
    nextChangeParam = 'b';
elseif strcmp(changeParam,'b')
    nextChangeParam = 'BGD';
end
end

