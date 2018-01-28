function eps=Lindhard(q,omega,gamma,omega0)
    
    n_dens = omega0*omega0 / (4.0*pi);
    E_f = 0.5*(3 * pi*pi*n_dens)^(2.0 / 3.0);
    v_f = (2 * E_f)^0.5;   % v_f=k_f in atomic units
     
    z = q / (2 * v_f);
    chi = sqrt(1.0 / (pi*v_f));
    
    z1 = omega / (q*v_f) + z +  1j*gamma / (q*v_f);
    d1 = vos_g(z1);
    z2 = omega / (q*v_f) - z +  1j*gamma / (q*v_f);
    d2 = vos_g(z2);
    
    epsreal = 1.0 + chi*chi / (4 * z*z*z)*(real(d1) - real(d2));
    epsimag = chi*chi / (4 * z*z*z)*(imag(d1) - imag(d2));
    eps = epsreal + epsimag*1j;

end