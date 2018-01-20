function eps=Mermin_LL(q,omega,gamma,omega0,omega_gap)
    
    g_over_w = gamma / omega;
    z1 = 1.0 + 1j*g_over_w; % omega should be unequal 0
    z2 = eps_LLX(q, omega, gamma, omega0, omega_gap) - 1.0;
    z3 = eps_LLX(q, 0.0, 0.0, omega0, omega_gap) - 1.0;
    
    top = z1*z2;
    bottom = 1.0 + 1j*g_over_w*z2 / z3;
    eps = 1.0 + top / bottom;
end