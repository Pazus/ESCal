function diimfp = nrm_diimfp_penn(ELF,E0)
%%
%{
   Calculates the normalised DIIMFP (in eV^-1) 
   from data for the dielectric loss function, i.e. 
   the imaginary part of the reciprocal of the dielectric function.
   The Penn algoritm is used.
   \param [in] w   - the energy loss array for which the diimfp is to be calculated
   \param [in] osc - the struct containing information on the oscilattor parameters
   \param [in] the Fermi energy
   \param [in] E0 - the energy for which the diimfp is to be calculated
%}
%%
    
    x = ELF(:,1);
    w = x(x<E0)/h2ev;
    
    x_in = zeros(length(w),1);
    
    for i = 1:length(w)
        if w(i)<0.05
            x_in(i) = 0;
        else
            hhw=w(i);
            wpmax=2*(hhw-E0/h2ev+sqrt(E0/h2ev*(E0/h2ev-hhw)));
            i1=0;
            k=2;
            quit=false;
            while ~quit
                if w(k)>wpmax
                    quit=true;
                    wup=wpmax;
                else
                    wup=w(k);
                end
                wlo=w(k-1);
                if wlo>wpmax
                    wlo=wpmax;
                end
                temp = ELF(k,2)*(ddmfpp_integrand(hhw,wup)-ddmfpp_integrand(hhw,wlo));
                if isnan(temp)
                    temp = 0;
                end
                i1 = i1 + temp;
                k = k+1;
                x_in(i) = i1/2.0/pi/hhw/a0/E0;
            end
        end
    end
    diimfp = x_in; % ./ trapz(osc.eloss,x_in);
end

function x = ddmfpp_integrand(de,hwp)
    if hwp>de
        error('ddmfpp_integrand:Somethings wrong');
    else
        x = -(hwp+de*log(de-hwp));
    end
end



