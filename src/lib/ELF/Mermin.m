function eps=Mermin(q, omega, gamma,omega0,alpha, isIonization)

    if nargin < 6
        isIonization=false;
    end
    
    sq = numel(q);
    q=q(:)';
    omega = omega(:);

    om_at_q = omega0 + 0.5*alpha*q.^2;
    
    g_over_w = gamma ./ omega;
    z1 = ones(size(g_over_w)) + 1j*g_over_w; % omega should be unequal 0
    z2 = Lindhard(q, omega, gamma, om_at_q) - complex(1,0);
    z3 = Lindhard(q, zeros(size(omega)), 0.0, om_at_q) - complex(1,0);
    
    top = bsxfun(@times,z1,z2);
    bottom = repmat(ones(size(omega)),1,sq) + bsxfun(@times,1j*g_over_w,z2)./z3;
    eps = repmat(ones(size(omega)),1,sq) + top./bottom;
    
    if isIonization
        ind = bsxfun(@lt,omega,om_at_q);
        eps(ind)=0;
    end
end