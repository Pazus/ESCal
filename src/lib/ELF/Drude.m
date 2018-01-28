function eps=Drude(q,w_global,omega0,gamma,alpha,FermiEnergy)

    if alpha < 1000
        % now there are different aplha's in use in the literature here we use dispersion relative to free particle
        %here not really important, dispersion is different than expected in Drude
        w_at_q = omega0 + 0.5*alpha * q*q;
    else
        % now we use full dispersion if alpha > 100 with the Fermi energy based on the total electron density, so no alpha fitting possible
        w_at_q_square = omega0 * omega0 + (2.0 / 3.0)*FermiEnergy*q*q + q*q*q *q/ 4.0;
        w_at_q = sqrt(w_at_q_square);
    end
    
    divisor = (w_global*w_global - w_at_q*w_at_q)*(w_global*w_global - w_at_q*w_at_q) + w_global*w_global*gamma * gamma;
	re_oneovereps =1.0+(omega0*omega0)*(w_global*w_global - w_at_q*w_at_q) / divisor;
	im_oneovereps = -(omega0*omega0)*w_global*gamma / divisor;
    
    eps = re_oneovereps + 1j*im_oneovereps;

end


