%% Function perturb_param

function PerturbParam = perturb_param(GenParam,noiseRatio)

pertR = max(diag(GenParam.R))*noiseRatio;
PerturbParam.R = GenParam.R + diag(pertR*randn(1,size(GenParam.R,2)));

[~,p] = chol(PerturbParam.R);
if p > 0
    while p > 0
        pertR = max(diag(GenParam.R))*noiseRatio;
        PerturbParam.R = GenParam.R + diag(pertR*randn(1,size(GenParam.R,2)));
        [~,p] = chol(PerturbParam.R);
    end
end

disp('R perturbed');
pertG = max(diag(GenParam.G))*noiseRatio;
PerturbParam.G = GenParam.G + diag(pertG*randn(1,size(GenParam.G,2)));

PerturbParam.G(PerturbParam.G>1) = 1;

pertB = max(diag(GenParam.B))*noiseRatio;
PerturbParam.B = GenParam.B - diag(pertB*rand(1,size(GenParam.B,2)));

pertQ = max(diag(GenParam.Q))*noiseRatio;
PerturbParam.Q = GenParam.Q + diag(pertQ*randn(1,size(GenParam.Q,2)));

[~,p] = chol(PerturbParam.Q);
if p > 0
    while p > 0
        pertQ = max(diag(GenParam.Q))*noiseRatio;
        PerturbParam.Q = GenParam.Q + diag(pertQ*randn(1,size(GenParam.Q,2)));
        [~,p] = chol(PerturbParam.Q);
    end
end
disp('Q perturbed');
pertD = max(diag(GenParam.D))*noiseRatio;
PerturbParam.D = GenParam.D + diag(pertD*randn(1,size(GenParam.D,2)));

pertA = max([GenParam.A(:)])*noiseRatio;
PerturbParam.A = GenParam.A + pertA*randn(size(GenParam.A,1),size(GenParam.A,2));

% pertV1 = max([GenParam.cov_1(:)])*noiseRatio;
% PerturbParam.cov_1 = GenParam.cov_1 + pertV1*randn(size(GenParam.cov_1,1),size(GenParam.cov_1,2));
% 
% [~,p] = chol(PerturbParam.cov_1);
% if p > 0
%     while p > 0
%         pertV1 = max([GenParam.cov_1(:)])*noiseRatio;
%         PerturbParam.cov_1 = GenParam.cov_1 + pertV1*randn(size(GenParam.cov_1,1),size(GenParam.cov_1,2));
%         [~,p] = chol(PerturbParam.cov_1);
%     end
% end
% disp('V_1 perturbed');
PerturbParam.cov_1 = GenParam.cov_1;

pertP = max(diag(GenParam.P))*noiseRatio;
PerturbParam.P = GenParam.P + diag(pertP*randn(1,size(GenParam.P,2)));

[~,p] = chol(PerturbParam.P);
if p > 0
    while p > 0
        pertP = max(diag(GenParam.P))*noiseRatio;
        PerturbParam.P = GenParam.P + diag(pertP*randn(1,size(GenParam.P,2)));
        [~,p] = chol(PerturbParam.P);
    end
end
disp('P perturbed');
pertd = max(GenParam.d)*noiseRatio;
PerturbParam.d = GenParam.d + pertd*randn(size(GenParam.d,2));

pertb = max(GenParam.b)*noiseRatio;
PerturbParam.b = GenParam.b + pertb*randn(size(GenParam.b,2));

pertk = max(GenParam.k)*noiseRatio;
PerturbParam.k = GenParam.k + pertk*randn(size(GenParam.k,2));

PerturbParam.mu_1 = zeros(size(GenParam.G,1),1);
end