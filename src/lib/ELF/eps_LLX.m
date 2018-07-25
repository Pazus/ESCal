function eps=eps_LLX(q,omega,gamma,omega0,omega_gap)
    
%     q=q(:)';
    omega = omega(:);
    
    ogdif = omega_gap^2 + gamma^2;

    omega_minus_square = complex(bsxfun(@minus,omega.^2,ogdif),2.0*omega*gamma);
    r = abs(omega_minus_square);
    theta = atan2(imag(omega_minus_square), real(omega_minus_square));
    omega_minus = complex(bsxfun(@times,sqrt(r),cos(theta/2.0)),bsxfun(@times,sqrt(r),sin(theta/2.0)));
    if bsxfun(@ge,real(omega_minus),0.0)
        eps = Lindhard(q,real(omega_minus),imag(omega_minus),omega0);
    else
        n_dens = omega0^2 / (4.0*pi);
        E_f = 0.5*(3.0 * pi^2*n_dens)^(2.0 / 3.0);
        v_f = (2 * E_f)^0.5;   % v_f=k_f in atomic units
        DeltaSquare = - omega_minus_square ./ (E_f^2);
        r = abs(DeltaSquare);
        theta = atan2(imag(DeltaSquare), real(DeltaSquare));
        Delta = bsxfun(@times,sqrt(r),cos(theta / 2.0)) + 1j*bsxfun(@times,sqrt(r),sin(theta / 2.0));
        QQ = q ./ v_f;
        z1 = 2.0*QQ+QQ.^2;
        res1 = z1 ./ Delta;
        res1 = c_arctan(res1);
        z2 = 2.0 * QQ + QQ.^2;
        res2 = bsxfun(@rdivide,z2, Delta);
        res2 = c_arctan(res2);
        res1 = res1 + res2;
        res2 = bsxfun(@times,res1,Delta);  % now res2= delta*(tan-1(2q+q^2/delta)+ tan-1(1(2q-q^2/delta)))
        
        z1 = bsxfun(@plus,real(DeltaSquare),bsxfun(@times,2 * QQ + QQ.^2,2 * QQ + QQ.^2)) + 1j*imag(DeltaSquare);
        z2 = bsxfun(@plus,real(DeltaSquare),bsxfun(@times,2 * QQ - QQ.^2,2 * QQ - QQ.^2)) + 1j*imag(DeltaSquare);
        z1 = z1 ./ z2;
        z1 = log(z1);
        z2 = bsxfun(@times,DeltaSquare,z1);
        
        p1 = bsxfun(@rdivide,imag(res2),2*QQ.^3);
        p2 = bsxfun(@rdivide,imag(z2),8*QQ.^5);
        p3 = bsxfun(@rdivide,imag(z1),2*QQ.^3);
        p4 = bsxfun(@rdivide,imag(z1),8*QQ);        
        eps_imag = bsxfun(@times,2.0 ./ (pi * v_f),-p1+p2+p3-p4);  
        
        t1 = bsxfun(@rdivide,real(res2),2*QQ.^3);
        t2 = bsxfun(@rdivide,real(z2),8*QQ.^5);
        t3 = bsxfun(@rdivide,real(z1),2*QQ.^3);
        t4 = bsxfun(@rdivide,real(z1),8*QQ);
        t5 = bsxfun(@minus,1.0./QQ.^2,t1);
        eps_real = ones(size(eps_imag)) + bsxfun(@times,2.0 ./ (pi * v_f),t5+t2+t3-t4);
        
        eps = complex(eps_real,eps_imag);
    end
end