function ELF = eps_sum(q,w,osc,FermiEnergy)

%%
%{
   \brief Calculates the matrix of the imaginary part of the inverse dielectric function
   for the energy w and momentum transfer q on the basis of four ELF models.
   \param [in] q - momentum transfer array (in 1/A)
   \param [in] w - energy loss array (in eV) for which Im(-1/eps) is to be calculated
   \param [in] osc - the struct containing information on the oscilattor parameters:
                            1. osc.A     - amplitudes (in eV^2)
                                sum(Ai) = 1 the systen is a mettalic;
                                sum(Ai) < 1 an insulator
                            2. osc.G     - damping coefficients (in eV) gamma
                            3. osc.Om    - plasmon (resonance) frequency omega
                            4. osc.alpha - alpha is a constant between 1 (metals) and 0 (insulators)
                            5. osc.u     - a quantity related to the band
                                gap (is used only for Mermin_LL)
                            6. osc.model - type of model: Drude, Lindhard, Mermin, Mermin_LL 
   \param [in] the Fermi energy
%}
%%


ELF = zeros(length(w),1);

    for k=1:length(q)
        for i=1:length(w)
            eps1 = 0;
            for j=1:length(osc.A)
                
                switch osc.model
                    case 'Drude'
                        epsDrud = Drude(q(k),w(i),osc.Om(j),osc.G(j),osc.alpha,FermiEnergy);
                        eps1 = eps1 + osc.A(j)*(epsDrud - 1.0);
                    case 'Lindhard'
                        epsMerm = Lindhard(q(k),w(i),osc.G(j),osc.Om(j));
                        eps1 = eps1 + osc.A(j)*(1.0/epsMerm - 1.0);
                    case 'Mermin'
                        epsMerm = Mermin(q(k),w(i),osc.G(j),osc.Om(j));
                        eps1 = eps1 + osc.A(j)*(1.0/epsMerm - 1.0);
                    case 'Mermin_LL'
                        epsMerm = Mermin_LL(q(k),w(i),osc.G(j),osc.Om(j),osc.u);
                        eps1 = eps1 + osc.A(j)*(1.0/epsMerm - 1.0);
                    otherwise
                        error('Choose the correct model name: Drude,Lindhard,Mermin,Mermin_LL');
                end                
            end  
            eps = 1.0/(eps1+1.0);
            ELF(i,k)= imag(-1/eps);
        end
    end    
end