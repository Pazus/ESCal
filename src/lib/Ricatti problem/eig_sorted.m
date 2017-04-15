function [V, D, invV, V_, invV_] = eig_sorted(X)
[V, D] = eig(X);
[D, iSort]=sort(diag(real(D)));
V = V(:,iSort);
if nargout>=3
    invV = inv(V);
end
if nargout==5
    V_ = V.';
    invV_ = invV.';
end
    