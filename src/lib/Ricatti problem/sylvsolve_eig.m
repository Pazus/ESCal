function X = sylvsolve_eig(a_eig, b_eig, C, U, V)
%SYLVSOLVE_EIG solve Lyapunov and Sylvester equations using
%Eigenvalues/eigenvectors factorization
%
%   A = U*a_eig*U^-1, B = V*b_eig*V^-1
%   a_eig, b_eig - eigenvalues of A and B correspondingly (diagonal matrices)
%   U,V - eigenvectors of A and B correspondingly
%
%   X = sylvsolve_eig(a_eig,C,U) - solves for Lyapunov equation
%
%        AX + XA' = C
%
%   X = sylvsolve_eig(a_eig,b_eig,C,U,V) - solves for Sylvester equation
%
%        AX + XB = C
%
%   Authors: Kaplya, Lubenchenko

    if nargin == 3 || (isequal(a_eig,b_eig) && isequal(U,V))
        U = C;
%         C = b_eig;
%         clear b_eig V;
%         C = U\C*U;
%         P = bsxfun(@plus,a_eig,a_eig');
%         X = U*(C./P)/U;
        C = U\b_eig/U';
        P = bsxfun(@plus,a_eig,a_eig');
        X = U*(C./P)*U';
    else
        C = U\C*V;
        P = bsxfun(@plus,a_eig,b_eig');
        X = U*(C./P)/V;
    end

end