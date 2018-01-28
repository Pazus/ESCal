function x_in=F_mod(delta,J_ion,lambda_ion,osc,E0,FermiEnergy)

delta_loc = delta;
x_pl = nrm_diimfp_penn(delta_loc,osc,FermiEnergy,E0);

x_ion_all = zeros(length(J_ion),length(delta_loc));

for i=1:length(J_ion)
    if ~(J_ion(i)>max(delta_loc))
        w_ion = 1 ./ (delta_loc.^2);
        w_ion(1:sum(delta_loc<J_ion(i))) = 0;
        x_ion = w_ion ./ trapz(delta_loc,w_ion);
        x_ion_all( i, : ) = x_ion;
    end
end
    
x_in_loc = x_pl';
for i=1:length(J_ion)
   x_in_loc = x_in_loc + lambda_ion(i)*(x_ion_all( i, : )-x_pl');
end

x_in = interp1(delta_loc,x_in_loc,delta);
end