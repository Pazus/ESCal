function eps=Mermin(q, omega, gamma, omega0)
    
    g_over_w = gamma / omega;
    z1 = 1.0 + g_over_w*1j; % omega should be unequal 0
    z2 = Lindhard(q, omega, gamma, omega0) - 1.0;
    z3 = Lindhard(q, 0.0, 0.0, omega0) - 1.0;
    
    top = z1*z2;
    bottom = 1.0 + g_over_w*1j*z2 / z3;
    eps = 1.0 + top / bottom;
end