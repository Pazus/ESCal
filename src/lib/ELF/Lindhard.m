function eps=Lindhard(q,omega,gamma,omega0)

    sq = numel(q);
    sw = numel(omega);
    
    q=q(:)';
    omega = omega(:);
    
    n_dens = bsxfun(@power,omega0,2)./ (4.0*pi);
    E_f = 0.5*(3 * pi^2*n_dens).^(2.0 / 3.0);
    v_f = bsxfun(@power,2*E_f,0.5);   % v_f=k_f in atomic units
     
    z = q./(2 * v_f);   
    chi = sqrt(1.0 ./ (pi*v_f));
    
    q_vf = bsxfun(@times,q,v_f);
    
    z1_1 = omega./q_vf;
    
    z1 = bsxfun(@plus,bsxfun(@plus,z1_1,z),1j*bsxfun(@rdivide,gamma,q_vf));
    d1 = vos_g(z1);
    z2 = bsxfun(@plus,bsxfun(@minus,z1_1,z),1j*bsxfun(@rdivide,gamma,q_vf));
    d2 = vos_g(z2);
    
    zzz = bsxfun(@power,z,3);
    red1_d2=bsxfun(@minus,real(d1),real(d2));
    imd1_d2=bsxfun(@minus,imag(d1),imag(d2));
    
    epsreal = repmat(ones(size(omega)),1,sq) + bsxfun(@times,red1_d2 ,bsxfun(@power,chi,2)./(4 * zzz));
    epsimag = bsxfun(@times,imd1_d2 ,bsxfun(@power,chi,2)./(4 * zzz));
    eps = epsreal + epsimag*1j;

end