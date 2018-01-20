function eps=eps_LLX(q,omega,gamma,omega0,omega_gap)

    omega_minus_square = omega*omega-omega_gap*omega_gap-gamma*gamma + 1j*2.0*omega*gamma;
    r = abs(omega_minus_square);
    theta = atan2(imag(omega_minus_square), real(omega_minus_square));
    omega_minus = sqrt(r)*cos(theta/2.0) + 1j*sqrt(r)*sin(theta/2.0);
    if (real(omega_minus) >= 0.0)
        eps = Lindhard(q,real(omega_minus),imag(omega_minus),omega0);
    else
        n_dens = omega0*omega0 / (4.0*pi);
        E_f = 0.5*(3.0 * pi*pi*n_dens)^(2.0 / 3.0);
        v_f = 2 * E_f^0.5;   % v_f=k_f in atomic units
        DeltaSquare = -omega_minus_square / (E_f*E_f);
        r = abs(DeltaSquare);
        theta = atan2(imag(DeltaSquare), real(DeltaSquare));
        Delta = sqrt(r)*cosd(theta / 2.0) + 1j*sqrt(r)*sind(theta / 2.0);
        QQ = q / v_f;
        z1 = 2.0*QQ+QQ*QQ + I*0;
        res1 = z1 / Delta;
        res1 = c_arctan(res1);
        z2 = 2.0 * QQ + QQ*QQ + I*0;
        res2 = z2 / Delta;
        res2 = c_arctan(res2);
        res1 = res1 + res2;
        res2 = res1*Delta;  % now res2= delta*(tan-1(2q+q^2/delta)+ tan-1(1(2q-q^2/delta)))
        
        z1 = (real(DeltaSquare) + (2 * QQ + QQ*QQ) * (2 * QQ + QQ*QQ)) + 1j*imag(DeltaSquare);
        z2 = (real(DeltaSquare) + (2 * QQ - QQ*QQ) * (2 * QQ - QQ*QQ)) + 1j*imag(DeltaSquare);
        z1 = z1 / z2;
        z1 = log(z1);
        z2 = DeltaSquare*z1;
        eps_imag =  2.0 / (pi * v_f)*(  - imag(res2) / (2*QQ*QQ*QQ)+imag(z2) /(8*QQ*QQ*QQ*QQ*QQ) + imag(z1) / (2 * QQ*QQ*QQ) - imag(z1) / (8 * QQ));
        eps_real = 1.0 + 2.0 / (pi * v_f)*(1.0/(QQ*QQ) - real(res2) / (2*QQ*QQ*QQ) +real(z2) /(8*QQ*QQ*QQ*QQ*QQ) + real(z1) / (2 * QQ*QQ*QQ) - real(z1) /(8 * QQ));
        eps = eps_real + 1j*eps_imag;
    end
end