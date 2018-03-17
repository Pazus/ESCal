function eps=Lindhard(q,omega,gamma,omega0)

    sq = numel(q);
    sw = numel(omega);
    
    q=q(:)';
    omega = omega(:);
    
    n_dens = omega0^2 / (4.0*pi);
    E_f = 0.5*(3 * pi^2*n_dens)^(2.0 / 3.0);
    v_f = (2*E_f)^0.5;   % v_f=k_f in atomic units
     
    z = q./(2 * v_f);   
    chi = sqrt(1.0 / (pi*v_f));
    
    z1_1 = omega./(q*v_f);
    
    z1 = bsxfun(@plus,bsxfun(@plus,z1_1,z),1j*gamma./(q*v_f));
    d1 = vos_g(z1);
    z2 = bsxfun(@plus,bsxfun(@minus,z1_1,z),1j*gamma./(q*v_f));
    d2 = vos_g(z2);
    
    zzz = z.^3;
    red1_d2 = real(d1) - real(d2);
    imd1_d2 = imag(d1) - imag(d2);
    
    epsreal = bsxfun(@plus,1,bsxfun(@times,red1_d2 ,chi^2./(4 * zzz)));
    epsimag = bsxfun(@times,imd1_d2 ,chi^2./(4 * zzz));
    eps = complex(epsreal,epsimag);

end