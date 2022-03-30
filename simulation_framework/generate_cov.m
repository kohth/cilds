%% Function generate_cov
%  Function generate_cov generates a symmetric positive definite matrix
%
%  Input:
%       dim - scalar indicating dimension of covariance matrix
%       noise - scalar indicating scale of covariance matrix
%
%  Output:
%       V - dimxdim symmetric positive definite matrix

function V = generate_cov(dim,noise)
%     V = noise*randn(dim);
%     V = (V*V')/2;
% %     [V,D,W] = eig(A);
%     [~,p] = chol(V);
%     if p > 0
%         warning("Matrix is not positive definite");
%     end
    
    Q = orth(randn(dim));

    lambda = rand(dim,1);
    lambda(lambda < .1) = .1;
    lambda = lambda*noise;
    V = Q*diag(lambda)*Q';
    [~,p] = chol(V);
    if p > 0
        warning('Matrix is not positive definite');
    end
end