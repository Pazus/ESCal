function R = care_Newton(A, B, C, D, eps,iter_max)
if nargin == 5
    iter_max = eps;
    eps = D;
    D = C;
    C = B;
    B = A.';
    symmetric = true;
else
    symmetric = (A == B');
end;
R = zeros(size(A));

R_old = ones(size(A));
iter = 1;
while iter < iter_max && norm(R-R_old)/norm(R)>eps
%     Gk = G-R*D;
    RD = R*D;
    if ~symmetric
        Ck = C - A*R - R*B + RD*R;
        R_old = R;
        R = R + lyap(A-RD,B-D*R,-Ck);
    else
        AR = A*R;
        ARD = A-RD;
        Ck = C - AR - AR' + RD*R;
        R_old = R;
        R = R + lyap(ARD,-Ck);
    end
    iter = iter + 1;
end
if iter >= iter_max
    error('No convergance!');
end
% disp(['Newton: ' num2str(iter) ' iterations']);
