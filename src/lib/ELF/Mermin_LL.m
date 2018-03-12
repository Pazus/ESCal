function eps=Mermin_LL(q,omega,gamma,omega0,omega_gap,alpha, isIonization)

    if nargin < 6
        isIonization=false;
    end

    sq = numel(q);
    sw = numel(omega);
    
    q=q(:)';
    omega = omega(:);
    
    om_at_q = omega0 + 0.5*alpha*q.^2;
    
    g_over_w = gamma ./ omega;
    z1 = ones(size(g_over_w)) + 1j*g_over_w; % omega should be unequal 0
    z2 = eps_LLX(q, omega, gamma, om_at_q, omega_gap) - ones(sw,sq);
    z3 = eps_LLX(q, zeros(size(omega)), 0.0, om_at_q, omega_gap) - ones(sw,sq);
    
    top = bsxfun(@times,z1,z2);
    bottom = ones(sw,sq) + bsxfun(@times,1j*g_over_w,z2)./z3;
    eps = ones(sw,sq) + top./bottom;
    
    if isIonization
        ind = bsxfun(@lt,omega,om_at_q);
        eps(ind)=0;
    end
end