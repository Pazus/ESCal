function diimfp = ndsep_Tung(osc,E0,theta,sign)
%%
%{
   Calculates the normalised NDSEP (in eV^-1) 
   for a given emission angle theta
   according to Tung algorithm
   C.J. Tung et al. / Phys. Rev. B 49 (1994) 16684.
   The Romberg algoritm for the integration is used.
%}
%%

x_in = zeros(length(osc.eloss),1);
w = osc.eloss;

for i=1:length(w) 
    qmin = sqrt(2*E0/h2ev)-sqrt(2*(E0/h2ev-w(i)/h2ev));
    qmax = sqrt(2*E0/h2ev)+sqrt(2*(E0/h2ev-w(i)/h2ev));
    sum = rombint_surf(osc,qmin/a0,qmax/a0,15,w(i),theta,E0,sign); %Romberg intergration
    x_in(i) = sum/(pi*E0);
    if ~mod(i,1000)
        YY = [num2str(i) ,' from ', num2str(length(w)-1)];
        disp(YY);
    end
end
diimfp = x_in ./ trapz(w,x_in);
%plot(osc.eloss,diimfp)
end