function ELF = eps_sum(osc)

%%
%{
   Calculates the imaginary part of the inverse dielectric function Im[-1/eps]
   for the energy w and momentum transfer q on the basis of four ELF models.
   Input parameter:
   osc - field with oscilattors data including:
         1. osc.A     - amplitudes
         2. osc.G     - damping coefficients gamma
         3. osc.Om    - plasmon (resonance) frequency omega
         4. osc.alpha - alpha is a constant between 1 (metals) and 0 (insulators)
         5. osc.u     - a quantity related to the band
                        gap (is used only for Mermin_LL)
         6. osc.model - type of model: Drude, Lindhard, Mermin, Mermin_LL
    In the case of Mermin and Mermin_LL models the array of q should always
    start from 0.01 (instead of 0).
%}
%%

if strcmp(osc.model,'Drude') 
    osc = convert2au(osc); %converts to atomic units, this is necessary for Drude model
end

elec_density = sum(osc.A)/4/pi;
Y = ['Electron density is ',num2str(elec_density),' (in a.u.)'];
disp(Y);

w = osc.eloss;
q = osc.qtran;

ELF = zeros(length(w),length(q));

    for k=1:length(q)
        for i=1:length(w)
            eps1 = complex(0.0,0.0);
            eps_re = osc.beps;
            eps_im = 0;
            for j=1:length(osc.A)
                
                switch osc.model
                    case 'Drude'
                        epsDrud_re = Drude_re(q(k),w(i),osc.Om(j),osc.G(j),osc.alpha,osc.Ef);
                        epsDrud_im = Drude_im(q(k),w(i),osc.Om(j),osc.G(j),osc.alpha,osc.Ef);
                        w_q = osc.Om(j) + 0.5*osc.alpha * q(k)*q(k);
                        if j>length(osc.A) - osc.numion && w(i)<w_q && osc.ion
                            eps_re = eps_re - 0;
                            eps_im = eps_im + 0;
                        else
                            eps_re = eps_re - osc.A(j)*epsDrud_re;
                            if (w(i)>osc.egap)
                                eps_im = eps_im + osc.A(j)*epsDrud_im;
                            else
                                eps_im = eps_im + 0;
                            end
                        end
                    case 'Lindhard'
                        w_at_q = osc.Om(j) + 0.5*osc.alpha * q(k)*q(k);
                        epsLind = Lindhard(q(k),w(i),osc.G(j),w_at_q);
                        eps1 = eps1 + osc.A(j)*(complex(1,0)/epsLind);
                    case 'Mermin'
                        w_at_q = osc.Om(j) + 0.5*osc.alpha * q(k)*q(k);
                        epsMerm = Mermin(q(k),w(i),osc.G(j),w_at_q);
                        eps1 = eps1 + osc.A(j)*(complex(1,0)/epsMerm);
                    case 'Mermin_LL'
                        w_at_q = osc.Om(j) + 0.5*osc.alpha * q(k)*q(k);
                        epsMerm = Mermin_LL(q(k),w(i),osc.G(j),w_at_q,osc.u);
                        eps1 = eps1 + osc.A(j)*(complex(1,0)/epsMerm);
                    otherwise
                        error('Choose the correct model name: Drude,Lindhard,Mermin,Mermin_LL');
                end                
            end 
            switch osc.model
                case 'Drude'
                    eps = complex(eps_re,eps_im);
                otherwise
                    eps = complex(1,0)/eps1;
            end
            ELF(i,k) = imag(-1/eps);
        end
    end 
end