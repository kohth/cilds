function [Xf,Vf,Vj, ll] = cifa_estep(Param, Observation)
% Perform expectation step of EM algorithm
% Refer to section 2.2 of supplemental methods for CIFA EM
%
% Example use:
%     [Xf, Vf, Vj, ll] = cifa_estep(Param,Observation)
%
% INPUT:
%
%      Param        - Structure containing model parameters
%                           -Dimensionality: 1 x 1
%                           -Fields:
%                               A (N_NEURONxN_LATENT),
%                               B (N_NEURONxN_NEURON),
%                               G (N_NEURONxN_NEURON) - Gamma,
%                               Q (N_NEURONxN_NEURON),
%                               R (N_NEURONxN_NEURON),
%                               P (N_NEURONxN_LATENT),
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
% OUTPUT:
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
%      ll           - Scalar containing computed loglikelihood of data
%
%% AUTHOR    : Koh Tze Hui
%% DEVELOPED : MATLAB (R2018a)
%% FILENAME  : cifa_estep.m
%% LAST CHECKED: 220331 (YYMMDD)
%% REFERENCES:
%     Yu B Derivation of Kalman Filtering and Smoothing Equations
%     Semedo J gLARA code
%%
compll = true;
speedup = false; %this looks for a convergence in K and J so it's an approximation
NProj = size(Param.A,1); NLatent = size(Param.A,2); 
R = Param.R;

% Reforming model
L = [Param.G, Param.A;zeros(NLatent,NProj),zeros(NLatent,NLatent)];
m = [Param.b; zeros(NLatent,1)];
S = [Param.Q, zeros(NProj,NLatent);zeros(NLatent,NProj),eye(NLatent,NLatent)];
nu_1 = [Param.mu_1; zeros(NLatent,1)];
sig_1 = [Param.cov_1, zeros(NProj,NLatent);zeros(NLatent,NProj),eye(NLatent,NLatent)];
Phi = [Param.B, zeros(NProj,NLatent)];

T = zeros(1,size(Observation,2)); NTrial = size(Observation,2);

Y = zeros(NTrial,NProj,size(Observation(1).y,2));
for iTrial = 1:NTrial
    T(iTrial) = size(Observation(iTrial).y,2);
    Y(iTrial,:,:) = Observation(iTrial).y;
end
ll = zeros(NTrial,1);
T = size(Observation(iTrial).y,2);
Xp = zeros(NTrial,NLatent+NProj,T);	% E[x_t|y_1,...,y_(t-1)]
Xc = zeros(NTrial,NLatent+NProj,T);	% E[x_t|y_1,...,y_t]
Xf = zeros(NTrial,NLatent+NProj,T);	% E[x_t|y_1,...,y_T]

Vp = zeros(NLatent+NProj,NLatent+NProj,T);	% cov[x_t|y_1,...,y_(t-1)]
Vc = zeros(NLatent+NProj,NLatent+NProj,T);	% cov[x_t|y_1,...,y_t]
Vf = zeros(NLatent+NProj,NLatent+NProj,T);	% cov[x_t|y_1,...,y_T]

K = zeros(NLatent+NProj,NProj,T);	% Kalman gain matrix
J = zeros(NLatent+NProj,NLatent+NProj,T);

Xp(:,:,1) = repmat(nu_1',[NTrial,1]) ;   % E[x_t]= mu_1, t=1
Vp(:,:,1) = sig_1;       % E[x_t*x_t^T] = V_1, t=1

invVy = inv(Phi*Vp(:,:,1)*Phi' + R); % Auxiliary computations
y_diff = squeeze(Y(:,:,1))-Xp(:,:,1)*Phi'; % Y-xPhi

Ktau        = Vp(:,:,1)*Phi'*invVy;
Xc(:,:,1)   = Xp(:,:,1) + y_diff*Ktau';
Vc(:,:,1)   = (eye(NLatent+NProj) - Ktau*Phi)*Vp(:,:,1);
Vc(:,:,1)   = (Vc(:,:,1) + Vc(:,:,1)')/2;

Xp(:,:,1+1) = Xc(:,:,1)*L'+[repmat(m',[NTrial,1])];
Vp(:,:,1+1) = L*Vc(:,:,1)*L'+S;
Vp(:,:,1+1) = (Vp(:,:,1+1) + Vp(:,:,1+1)')/2;

if compll
    % Update log-likelihood
    logdetinvS = logdet(invVy);
    ll = ll + (NTrial/2)*logdetinvS - 0.5*sum( sum( y_diff.*(y_diff*invVy) ) );
end


invR    = inv(R);
invR    = (invR + invR')/2;
PhiinvR   = Phi'*invR;

%% Kalman filter (forward step)
for t = 1+1:T
    y_diff = Y(:,:,t)-Xp(:,:,t)*Phi'; 
    if speedup && t > 2 && sum(sum((squeeze(K(:,:,t-1))-squeeze(K(:,:,t-2))<1e-15)))==size(K,1)*size(K,2)
        K(:,:,t) = K(:,:,t-1);
    else
        K(:,:,t)    = Vp(:,:,t)*Phi'/(R+Phi*Vp(:,:,t)*Phi');
        invVy = invR*(eye(NProj) - Phi*K(:,:,t));
        invVy = (invVy + invVy')/2;
    end
    Xc(:,:,t)   = Xp(:,:,t) + y_diff*K(:,:,t)';
    Vc(:,:,t)   = (eye(NLatent+NProj) - K(:,:,t)*Phi)*Vp(:,:,t);
    Vc(:,:,t)   = (Vc(:,:,t) + Vc(:,:,t)')/2;
    
    if t < T
        Xp(:,:,t+1) = Xc(:,:,t)*L'+[repmat(m',[NTrial,1])];
        Vp(:,:,t+1) = L*Vc(:,:,t)*L'+S;
        Vp(:,:,t+1) = (Vp(:,:,t+1) + Vp(:,:,t+1)')/2;
    end
    if compll
        % Update log-likelihood
        logdetinvS = logdet(invVy);
        ll = ll + (NTrial/2)*logdetinvS - 0.5*sum( sum( y_diff.*(y_diff*invVy) ) );
    end
end
if compll; ll = ll - 0.5*NTrial*NProj*T*log(2*pi); end

%% Kalman smoother (backwards step)
Xf(:,:,T) = Xc(:,:,T);
Vf(:,:,T) = Vc(:,:,T);
for t = 1:T-1
    if speedup && t>3 && sum(sum(abs(squeeze(J(:,:,t-1))-squeeze(J(:,:,t-2)))<1e-14))==size(J,1)*size(J,2)
        J(:,:,t) = J(:,:,t-1);
    else
        J(:,:,t) = Vc(:,:,t)*L'/(Vp(:,:,t+1));
    end
end
for t = (T-1):-1:1
    Xf(:,:,t) = Xc(:,:,t) + (Xf(:,:,t+1)-Xc(:,:,t)*L'-[repmat(m',[NTrial,1])])*J(:,:,t)';
    Vf(:,:,t) = Vc(:,:,t) + J(:,:,t)*(Vf(:,:,t+1) - Vp(:,:,t+1))*J(:,:,t)';
    Vf(:,:,t) = (Vf(:,:,t) + Vf(:,:,t)')/2;
    Vj(:,:,t+1) = Vf(:,:,t+1)*J(:,:,t)';
end
ll = mean(ll);
end




