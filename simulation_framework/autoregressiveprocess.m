function z = autoregressiveprocess(Param,tau,varargin)

p = inputParser;
p.addRequired('Param',@isstruct);
p.addRequired('tau',@(x) isnumeric(x));
p.parse(Param,tau,varargin{:});

Param = p.Results.Param;
tau = p.Results.tau;
tauInd = tau(1)==Param.TAU;
T = Param.T_GEN; p = Param.N_LATENT;

z(:,1) = randn(p,1);
D = converttautod(tau);
% D = Param.DCONST(tauInd)*ones(p,1);
D = diag(D);
P = eye(p)-D*D';
for i = 2:T
    z(:,i) = D*z(:,i-1)+mvnrnd(zeros(p,1),P)';
end

% figure();
% subplot(2,1,1)
% plot(z(1,:),'k','LineWidth',1.5);
% subplot(2,1,2)
% plot(z(2,:),'k','LineWidth',1.5);