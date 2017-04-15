function X = symmetrize(X)
s = size(X);
if s(1) ~= s(2)
    error('Matrix has to be square');
end

if size(X,3)>1
    for i=1:size(X,3)
        X(:,:,i) = (X(:,:,i) + X(:,:,i)')/2;
    end
else
    X = (X + X')/2;
end

% if nargin == 1 || (upper ~= 1 && upper ~= false)
%     idxL = logical( tril(ones(size(X)),-1) );
%     X2 = X.';
%     X(idxL) = X2(idxL);
% else
%     idxL = logical( triu(ones(size(X)),1) );
%     X2 = X.';
%     X(idxL) = X2(idxL);
% end