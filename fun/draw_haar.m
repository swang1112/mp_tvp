function out = draw_haar(K,n)
%draw_haar draw orthonormal matrix from Haar prior
%   K: dimension of space (row of Q)
if n == K
    X = randn(K, K);
    [Q,R] = qr(X);
    Q(:,diag(R)<0) = Q(:,diag(R)<0)*-1;
    out = Q;
elseif n == 1
    X = randn(K, 1);
    out = X / norm(X);
end
        
end

