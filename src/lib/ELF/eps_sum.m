function ELF = eps_sum(q,w,osc,FermiEnergy)
%Perform the 


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