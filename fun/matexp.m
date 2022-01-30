function out = matexp(X, n)
%matexp matrix exponentials
if n == 0
    out = eye(size(X));
elseif n == 1
    out = X;
else
    out = X*matexp(X, n-1);
end
end
