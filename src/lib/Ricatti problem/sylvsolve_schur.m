% solving AX + BX = C;
% U*A*U', V*B*V' - Schur factorization
function X = sylvsolve_schur(R1,R2,C0,U,V)
    m = size(R1,2);
    n = size(R2,1);
    C = U'*C0*V;
    Is = eye(m);
    X = zeros(m,n);

    if any(diag(R2,-1)~=0)
        if nargin < 6
            Usq = U*U;
        end;
        i=1;
        while i < n + 1
            if i < n && R2(i+1,i)>10^-12*max(abs(R2(i,i)),abs(R2(i+1, i+1)))

                r11 = R2(i,i); r12 = R2(i,i+1);
                r21 = R2(i+1,i); r22 = R2(i+1,i+1);

                b = C(:,i:i+1)-X(:,1:i-1)*R2(1:i-1, i:i+1);
                b = [R1*b(:,1)+r22*b(:,1)-r21*b(:,2), R1*b(:,2)+r11*b(:,2)-r12*b(:,1)];
                X(:,i:i+1) = (Usq+(r11+r22)*R1+(r11*r22-r12*r21)*Is)\b;
                i = i+2;
            else
                b = C(:,i)-X(:,1:i-1)*R2(1:i-1,i);
                X(:,i) = (R1+R2(i,i)*Is)\b;
                i = i+1;
            end
        end
    else
        for jj=1:n
            b = C(:,jj)-X(:,1:jj-1)*R2(1:jj-1,jj);
            X(:,jj) = (R1+R2(jj,jj)*Is)\b;
        end
    end;

    X = U*X*V';
end