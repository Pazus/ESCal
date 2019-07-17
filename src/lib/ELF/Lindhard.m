function [eps,epsimag]=Lindhard(q,omega,gamma,omega0)

%     sq = numel(q);
%     sw = numel(omega);
    
%     q=q(:)';
    omega = omega(:);
    
    n_dens = omega0^2 / (4.0*pi);
    E_f = 0.5*(3 * pi^2*n_dens)^(2.0 / 3.0);
    v_f = (2*E_f)^0.5;   % v_f=k_f in atomic units
     
    z = q./(2 * v_f);   
    chi = sqrt(1.0 / (pi*v_f));
    
    z1_1 = omega./(q*v_f);
    z1_1(isnan(z1_1)) = 0;
    
    gq = gamma./(q*v_f); gq(isnan(gq)) = 0;
    [reD1, imD1] = vos_g(z1_1 + z, gq);
    [reD2, imD2] = vos_g(z1_1 - z, gq);
    
    red1_d2 = reD1 - reD2;
    imd1_d2 = imD1 - imD2;
    
    chizzz = chi^2./(z.*z.*z * 4);
    epsreal = 1 + red1_d2 .* chizzz;
    epsimag = imd1_d2 .* chizzz;
    if nargout == 1
        eps = complex(epsreal,epsimag);
    else
        eps = epsreal;
        % epsimag = epsimag
    end

end