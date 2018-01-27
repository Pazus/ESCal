function ELF = ELF_k(q,w,osc,FermiEnergy)

%%
%{
   \brief Calculates the imaginary part of the inverse dielectric function
   for the energy loss array w and momentum transfer q on the basis of four ELF models.
   \param [in] q momentum transfer (in 1/A)
   \param [in] w energy loss array (in eV)
   \param [in] osc field with oscilattors data including:
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
elec_density = 0;
for i = 1:length(osc.A)
    elec_density = elec_density + osc.A(i) * osc.Om(i) * osc.Om(i)/(4.0 * pi); %Bethe sum rule
end
m = 0.511e6; %eV/(c^2)
h = 1240; %eV*nm
N = 8*sqrt(2)*pi*m^(3/2)*(2/3*FermiEnergy^(3/2))/(h^3); %the free electron density of a metal
X = ['Bethe sum rule result: elec_density is ',num2str(elec_density)];
disp(X)
Y = ['The free electron density is ',num2str(N),'/nm^3'];
disp(Y)

    for k=1:length(w)
            eps1 = 0;
            for j=1:length(osc.A)
                
                switch osc.model
                    case 'Drude'
                        epsDrud = Drude(q,w(k),osc.Om(j),osc.G(j),osc.alpha,FermiEnergy);
                        eps1 = eps1 + osc.A(j)*(epsDrud - 1.0);
                    case 'Lindhard'
                        epsMerm = Lindhard(q,w(k),osc.G(j),osc.Om(j));
                        eps1 = eps1 + osc.A(j)*(1.0/epsMerm - 1.0);
                    case 'Mermin'
                        epsMerm = Mermin(q,w(k),osc.G(j),osc.Om(j));
                        eps1 = eps1 + osc.A(j)*(1.0/epsMerm - 1.0);
                    case 'Mermin_LL'
                        epsMerm = Mermin_LL(q,w(k),osc.G(j),osc.Om(j),osc.u);
                        eps1 = eps1 + osc.A(j)*(1.0/epsMerm - 1.0);
                    otherwise
                        error('Choose the correct model name: Drude,Lindhard,Mermin,Mermin_LL');
                end                
            end  
            eps = 1.0/(eps1+1.0);
            ELF(k)= imag(-1/eps);
    end    
end