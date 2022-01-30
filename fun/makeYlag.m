function out = makeYlag(y,p)
%makeYlag make the matrix of lagged y
%   y is of dimension Tx1, the output is of dimension (T-p)xp
nrows = size(y,1);
YLag = zeros(nrows, p);
for i = 1:p
    YLag(i:nrows, i) = y(1:(nrows-i+1));
end
out = YLag(p:(nrows-1), 1:p);
end

