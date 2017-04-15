function R = care_MSS(A, D, C0, K)
%$MSS - solves CARE using MSS method
%
%       A*R + R*A' = C0 + R*D*R
%
%   K - number of strong scatterings

[U, a_eig] = eig_sorted(A);
[V, b_eig] = eig_sorted(A');
s = size(A);
R_temp = zeros(s(1),s(2),K+1);
% R_temp(:,:,1) = sylvsolve_eig(a_eig,C0,U);
R_temp(:,:,1) = sylvsolve_eig(a_eig,b_eig,C0,U,V);

for i=2:K+1
    C = zeros(s);
    for j=1:floor(i/2)
        C = C + R_temp(:,:,j)*D*R_temp(:,:,i-j);
    end
    C = C + C';
    if iseven(i)
        C = C + R_temp(:,:,i/2)*D*R_temp(:,:,i/2);
    end
    R_temp(:,:,i) = sylvsolve_eig(a_eig,C,U);
end

R = sum(R_temp,3);
