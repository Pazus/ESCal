function eps=Mermin_LL(q,omega,gamma,omega0,omega_gap)
    
%     q=q(:)';
    omega = omega(:);
    
    om_at_q = omega0;
    
    g_over_w = gamma ./ omega;
    z1 = complex(1,g_over_w); % omega should be unequal 0
    z2 = bsxfun(@minus,eps_LLX(q, omega, gamma, om_at_q, omega_gap),1);
    z3 = bsxfun(@minus,eps_LLX(q, zeros(size(omega)), 0.0, om_at_q, omega_gap),1);
    
    top = bsxfun(@times,z1,z2);
    bottom = bsxfun(@plus,1,bsxfun(@times,complex(0,g_over_w),z2)./z3);
    eps = bsxfun(@plus,1,top./bottom);
    
end