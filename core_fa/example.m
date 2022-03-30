% Generate fake data using FA model
% Data dimensionality: 20
% Latent dimensionality: 5
rng 'default'
params.L  = randn(20, 1);
params.Ph = 0.1 * (1:20)';
params.mu = ones(20, 1);

[Z,X] = simdata_fa(params, 100);

% Cross-validation
dim = crossvalidate_fa(X);

% % Identify optimal latent dimensionality
istar = ([dim.sumLL] == max([dim.sumLL]));

% Project training data into low-d space
[estParams, LL] = fastfa(X, 1);

Z_test = fastfa_estep(X, estParams);

figure();
plot(Z');

figure();
plot(-Z_test.mean');