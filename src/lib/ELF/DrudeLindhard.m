function [eps_re, eps_im]=DrudeLindhard(q,w_global,omega0,gamma,alpha,FermiEnergy)

    w_global = w_global(:);

    if alpha < 1000
        % now there are different aplha's in use in the literature here we use dispersion relative to free particle
        %here not really important, dispersion is different than expected in Drude
        %w_at_q = omega0 + 0.5*alpha * q*q;
        w_at_q = omega0 + 0.5*alpha * q.^2;
    else
        % now we use full dispersion if alpha > 100 with the Fermi energy based on the total electron density, so no alpha fitting possible
        %w_at_q_square = omega0 * omega0 + (2.0 / 3.0)*FermiEnergy*q*q + q*q*q *q/ 4.0;
        w_at_q_square = omega0 * omega0 + (2.0 / 3.0)*FermiEnergy*q.^2 + q.^4/ 4.0;
        w_at_q = sqrt(w_at_q_square);
    end
    
    mm = bsxfun(@minus,w_global.^2,w_at_q.^2);
    divisor = mm.^2 + w_global.^2*gamma^2;

    eps_re = bsxfun(@minus,w_global.^2,w_at_q.^2)*(omega0.^2)./divisor;
    eps_im = gamma*(omega0.^2).*bsxfun(@rdivide,w_global,divisor);
    
end