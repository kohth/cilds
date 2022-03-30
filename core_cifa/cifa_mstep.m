function [EstParam,nextChangeParam] = cifa_mstep(Xf,Vf,Vj,PrevParam,Observation,FixParam,changeParam)
% Perform maximization step of EM algorithm
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
%                               Q (N_NEURONxN_NEURON),
%                               R (N_NEURONxN_NEURON),
%                               b (N_NEURONx1),
%                               mu_1 (N_NEURONx1),
%                               cov_1 (N_NEURONxN_NEURON)
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
%                               A,B,G,Q,R,b,mu_1,cov_1
% 
%      changeParam - String indicating which parameters to maximize in this
%                    partial maximization step
%
% OUTPUT:
%      EstParam    - Structure containing model parameter after current
%                    maximization step
%                           - Dimensions: 1 x 1
%                           - Fields:
%                               A,B,G,Q,R,b,mu_1,cov_1
%
%      nextChangeParam - String indicating parameters to maximize in next
%                        partial maximization step
%                           - Groupings: 'BG','A','b'
%% AUTHOR    : Koh Tze Hui
%% DEVELOPED : MATLAB (R2018a)
%% FILENAME  : cifa_mstep.m
%% LAST CHECKED: 220331 (YYMMDD)
%% REFERENCES:
%     Ghahramani Z Parameter Estimation for Linear Dynamical Systems


NLatent = size(PrevParam.A,2); NProj = size(Observation(1).y,1);
NTrial = size(Observation,2);
T = size(Observation(1).y,2);


%% Initialization for summations
sumyc = zeros(NProj,NProj);
sumyy = zeros(NProj,NProj);
sumcc_1 = zeros(NProj,NProj);
sumzc_1 = zeros(NLatent,NProj);
sumcz = zeros(NProj,NLatent);
sumz1z1 = zeros(NLatent,NLatent);
sumcc = zeros(NProj,NProj);
sumz1z = zeros(NLatent,NLatent);
cV = Vf(1:NProj,1:NProj,1:end);
zV = Vf(NProj+1:NProj+NLatent,NProj+1:NProj+NLatent,:);
cc_1 = Vj(1:NProj,1:NProj,:);
cz1 = Vf(1:NProj,NProj+1:NProj+NLatent,:);
cz = Vj(1:NProj,NProj+1:NProj+NLatent,:);
z1z = Vj(NProj+1:NProj+NLatent,NProj+1:NProj+NLatent,:);
cfut = Xf(:,1:NProj,:);
z1fut = Xf(:,NProj+1:NProj+NLatent,:);

y = zeros(NTrial,NProj,size(Observation(1).y,2));
for iTrial = 1:NTrial
    y(iTrial,:,:) = Observation(iTrial).y;
end
%% Summations

for t = T(1):-1:1
    sumcc =sumcc + cV(:,:,t) + cfut(:,:,t)'*cfut(:,:,t)/NTrial;
    sumyc = sumyc + y(:,:,t)'*cfut(:,:,t)/NTrial;
    sumyy = sumyy + y(:,:,t)'*y(:,:,t)/NTrial;
    sumz1z1 = sumz1z1+ zV(:,:,t) + z1fut(:,:,t)'*z1fut(:,:,t)/NTrial;
    if t > 1
        sumcc_1 = sumcc_1 + cc_1(:,:,t) + cfut(:,:,t)'*cfut(:,:,t-1)/NTrial;
        sumzc_1 = sumzc_1 + cz1(:,:,t-1)' + z1fut(:,:,t-1)'*cfut(:,:,t-1)/NTrial;
        sumcz = sumcz + cz(:,:,t) + cfut(:,:,t)'*z1fut(:,:,t-1)/NTrial;
        sumz1z = sumz1z  +z1z(:,:,t) + z1fut(:,:,t)'*z1fut(:,:,t-1)/NTrial;
    end
end

sumc_1c_1 = sumcc - cV(:,:,T) - cfut(:,:,T)'*cfut(:,:,T)/NTrial;
sumzz = sumz1z1 - zV(:,:,T) - z1fut(:,:,T)'*z1fut(:,:,T)/NTrial;
sumz1z1 = sumz1z1 - zV(:,:,1) - z1fut(:,:,1)'*z1fut(:,:,1)/NTrial;

sumc_1c_1 = 0.5*(sumc_1c_1+sumc_1c_1');
sumzz = 0.5*(sumzz+sumzz');
sumcc = 0.5*(sumcc+sumcc');


sumc = sum(sum(cfut,1)/NTrial,3)';
sumc_1 = sum(sum(cfut(:,:,1:end-1),1)/NTrial,3)';
sumz = sum(sum(z1fut(:,:,1:end-1),1)/NTrial,3)';
%% B, d, R

if isfield(FixParam,'B')
    B = FixParam.B;
else
    B = PrevParam.B;
end


if ~isfield(FixParam,'B')
    if strcmp(changeParam,'BG')
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
    if strcmp(changeParam,'BG')
        G = (diag(diag(sumcc_1)) - diag(diag(A*sumzc_1)) - diag(diag(b*sumc_1')))/(diag(diag(sumc_1c_1)));
    end
end

if ~isfield(FixParam,'A')
    if strcmp(changeParam,'Adk')
        A = (sumcz - G*sumzc_1' - b*sumz')/(sumzz);
    end
end

if ~isfield(FixParam,'b')
    if strcmp(changeParam,'b')
        b = (1/(T-1))*(sumc -(sum(cfut(:,:,1),1)/NTrial)' - G*sumc_1 - A*sumz);
    end
end

EstParam.b = b;
EstParam.A = A;
EstParam.G = G;

if ~isfield(FixParam,'Q')
    Q = (1/(T-1))*(diag(diag(sumcc - cV(:,:,1) - cfut(:,:,1)'*cfut(:,:,1)/NTrial))-2*diag(diag(G*sumcc_1'))...
        -2*diag(diag(A*sumcz'))-2*diag(diag(b*(sumc-(sum(cfut(:,:,1),1)/NTrial)')'))+2*diag(diag(A*sumzc_1*G'))...
        +diag(diag(G*sumc_1c_1*G'))+2*diag(diag(b*sumc_1'*G'))+diag(diag(A*sumzz*A'))...
        +2*diag(diag(b*sumz'*A'))+(T-1)*diag(diag(b*b')));
else
    Q = FixParam.Q;
end
EstParam.Q = Q;

%% mu_1, cov_1
if ~isfield(FixParam,'mu_1')
    mu_1 = (sum(cfut(:,:,1),1))';
    mu_1 = mu_1/NTrial;
else
    mu_1 = FixParam.mu_1;
end

if ~isfield(FixParam,'cov_1')
    cov_1 = diag(diag(cV(:,:,1))) + diag(diag((cfut(:,:,1)-repmat(mu_1,1,NTrial)')'*(cfut(:,:,1)-repmat(mu_1,1,NTrial)')))/NTrial;
else
    cov_1 = FixParam.cov_1;
end
cov_1 = 0.5*(cov_1+cov_1');
EstParam.cov_1 = cov_1;
EstParam.mu_1 = mu_1;

if strcmp(changeParam,'BG')
    nextChangeParam = 'A';
elseif strcmp(changeParam,'A')
    nextChangeParam = 'b';
elseif strcmp(changeParam,'b')
    nextChangeParam = 'BG';
end
end

