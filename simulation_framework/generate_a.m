%% Function generate_a
%  Function generate_a creates a square matrix (dim PxP where P is the 
%  dimension of the latent) where the eigenvalues are less than 1 to avoid
%  unstable systems which blow up
%
%  Input:
%       dim - scalar indicating the latent (x) dimension
%       coeff - scalar to scale A
%
%  Output:
%       A - a dimxdim matrix with eigenvalues less than 1

function A =generate_a(nProj,nLatent,coeff)
if nargin <1
    coeff = 10;
    nProj = 100;
    nLatent = 10;
    
end
    Q = orth(randn(nProj));
    W = orth(randn(nLatent));

    lambda = zeros(nProj,nLatent);
%     x = linspace(0.1,1,nLatent);
%     k = 5.0;
%     f = @(x,k) exp(-abs(k*x));
    load('V1_faparams.mat');
    eigspec = svd(faParams.L);
    lambda(1:nLatent,1:nLatent) = diag(coeff*eigspec);
    A = Q*lambda*W';

end