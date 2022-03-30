function [Xf,Vf,Vj, ll] = lds_estep(Param, Observation)
% Perform expectation step of EM algorithm
%
% Example use:
%     [Xf, Vf, Vj, ll] = lds_estep(Param,Observation)
%
% INPUT:
%
%      Param        - Structure containing model parameters
%                           -Dimensionality: 1 x 1
%                           -Fields:
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
% OUTPUT:
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
%      ll           - Scalar containing computed loglikelihood of data
%
%% AUTHOR    : Koh Tze Hui
%% DEVELOPED : MATLAB (R2018a)
%% FILENAME  : lds_estep.m
%% LAST CHECKED: 211101 (YYMMDD)
%% REFERENCES:
%     Yu B Derivation of Kalman Filtering and Smoothing Equations
%     Semedo J gLARA code
%% === Assign model parameters ===
compll = true; 
D = Param.D; P = Param.P; R = Param.R;
b = Param.b; A = Param.A; 
h_1 = Param.h_1; G_1 = Param.G_1;

%% === Reformat structure into matrix ===
T = zeros(1,size(Observation,2)); N_TRIAL = size(Observation,2);
for iTrial = 1:N_TRIAL
    T(iTrial) = size(Observation(iTrial).y,2);
    Y(iTrial,:,:) = Observation(iTrial).y;
end
%% === Create variables ===
N_NEURON = size(A,1); N_LATENT = size(A,2);

ll = zeros(N_TRIAL,1);
T = size(Observation(iTrial).y,2);
Xp = zeros(N_TRIAL,N_LATENT,T);	% E[z_t|y_1,...,y_(t-1)]
Xc = zeros(N_TRIAL,N_LATENT,T);	% E[z_t|y_1,...,y_t]
Xf = zeros(N_TRIAL,N_LATENT,T);	% E[z_t|y_1,...,y_T]

Vp = zeros(N_LATENT,N_LATENT,T);	% cov[z_t|y_1,...,y_(t-1)]
Vc = zeros(N_LATENT,N_LATENT,T);	% cov[z_t|y_1,...,y_t]
Vf = zeros(N_LATENT,N_LATENT,T);	% cov[_t|y_1,...,y_T]

K = zeros(N_LATENT,N_NEURON,T);	% Kalman gain matrix
J = zeros(N_LATENT,N_LATENT,T);

%% === Initialize variables at t = 1 ===
Xp(:,:,1) = repmat(h_1',[N_TRIAL,1]); % E[z_t]= h_1, t=1
Vp(:,:,1) = G_1;  % E[z_t*z_t^T] = G_1, t=1

invVy = inv(A*Vp(:,:,1)*A' + R); % Auxiliary computations
y_diff = (squeeze(Y(:,:,1))-[ repmat(b',[N_TRIAL,1])] )-...
    Xp(:,:,1)*A'; % Y-b-xA
Ktau        = Vp(:,:,1)*A'*invVy;
Xc(:,:,1)   = Xp(:,:,1) + y_diff*Ktau';
Vc(:,:,1)   = (eye(N_LATENT) - Ktau*A)*Vp(:,:,1);
Vc(:,:,1)   = (Vc(:,:,1) + Vc(:,:,1)')/2;                             

Xp(:,:,1+1) = Xc(:,:,1)*D';
Vp(:,:,1+1) = D*Vc(:,:,1)*D'+P;
Vp(:,:,1+1) = (Vp(:,:,1+1) + Vp(:,:,1+1)')/2;    

invR    = inv(R);
invR    = (invR + invR')/2;
AinvR   = A'*invR;
AinvRA  = AinvR*A;

%% === Compute log-likelihood ===
if compll
    logdetinvS = logdet(invVy);
    ll = ll + (N_TRIAL/2)*logdetinvS - 0.5*sum( sum( y_diff.*(y_diff*invVy) ) );
end

%% === Kalman filter (forward step)===
for t = 1+1:T
    y_diff = (Y(:,:,t)-repmat(b',[N_TRIAL,1]))-Xp(:,:,t)*A';      % Auxiliary computations    
    K(:,:,t)	= ((eye(N_LATENT) + Vp(:,:,t)*AinvRA)\Vp(:,:,t))*AinvR;
    Xc(:,:,t)   = Xp(:,:,t) + y_diff*K(:,:,t)';
    Vc(:,:,t)   = (eye(N_LATENT) - K(:,:,t)*A)*Vp(:,:,t);
    Vc(:,:,t)   = (Vc(:,:,t) + Vc(:,:,t)')/2;                               
        
    invVy = invR*(eye(N_NEURON) - A*K(:,:,t)); % Auxiliary computations
    invVy = (invVy + invVy')/2;
    
    if t < T
        Xp(:,:,t+1) = Xc(:,:,t)*D';
        Vp(:,:,t+1) = D*Vc(:,:,t)*D'+P;
        Vp(:,:,t+1) = (Vp(:,:,t+1) + Vp(:,:,t+1)')/2;                       
    end

    if compll
        % Update log-likelihood
        logdetinvS = logdet(invVy);
        ll = ll + (N_TRIAL/2)*logdetinvS - 0.5*sum( sum( y_diff.*(y_diff*invVy) ) );
    end
end
if compll; ll = ll - 0.5*N_TRIAL*N_NEURON*T*log(2*pi); end

%% === Kalman smoother (backwards step) ===
Xf(:,:,T) = Xc(:,:,T);
Vf(:,:,T) = Vc(:,:,T);
for t = (T-1):-1:1
    J(:,:,t) = Vc(:,:,t)*D'/(Vp(:,:,t+1));
    
    Xf(:,:,t) = Xc(:,:,t) + (Xf(:,:,t+1)-Xc(:,:,t)*D')*J(:,:,t)';
    Vf(:,:,t) = Vc(:,:,t) + J(:,:,t)*(Vf(:,:,t+1) - Vp(:,:,t+1))*J(:,:,t)';
    Vf(:,:,t) = (Vf(:,:,t) + Vf(:,:,t)')/2;    
    Vj(:,:,t+1) = Vf(:,:,t+1)*J(:,:,t)';
end


ll = mean(ll);
Xf = [Xf,ones(N_TRIAL,1,T)]; % Add a row for the constant in m-step
Vf = [[Vf,zeros(N_LATENT,1,T)];zeros(1,N_LATENT+1,T)];
Vj = [[Vj,zeros(N_LATENT,1,T)];zeros(1,N_LATENT+1,T)];






