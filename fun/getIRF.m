function out = getIRF(A,B,h)
%getIRF get IRFs
%   A: AR parameters, B: structural impact multipliers, h: horizons

[K, Kp] = size(A);
p = Kp/K;
out = zeros(K, K, h);
out(1) = B;
if p == 1
    for i = 2:h
        out[i] = matexp(A, i-1) * B;
    end
else
    Mm = zeros(Kp, Kp);
    Mm(1:K, 1:Kp) = A;
    Mm((K+1):Kp, 1:(p - 1) * K) = eye(K * (p - 1), K * (p - 1));
    Mm1 = eye(Kp, Kp);
    for i = 2:h
        Mm1 = Mm1 * Mm;
        out[i] = Mm1(1:K, 1:K)*B;
    end
end
    
   
