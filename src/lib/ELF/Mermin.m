function eps=Mermin(q, omega, gamma,omega0, isIonization)

    if nargin < 5
        isIonization=false;
    end
    
    sq = numel(q);
    sw = numel(omega);
    
    q=q(:)';
    omega = omega(:);
    
    om_at_q = omega0;
    
    g_over_w = gamma ./ omega;
    z1 = ones(size(g_over_w)) + 1j*g_over_w; % omega should be unequal 0
    z2 = Lindhard(q, omega, gamma, om_at_q) - ones(sw,sq);
    z3 = Lindhard(q, zeros(size(omega)), 0.0, om_at_q) - ones(sw,sq);
    
    top = bsxfun(@times,z1,z2);
    bottom = ones(sw,sq) + bsxfun(@times,1j*g_over_w,z2)./z3;
    eps = ones(sw,sq) + top./bottom;
    
%     if isIonization
%         for i=1:sq
%             [omax posmax]=max(imag(-1./eps(:,i)));
%             ind = bsxfun(@lt,omega,omega(posmax));
%             eps(ind,i)=0;
%         end
%     end
end