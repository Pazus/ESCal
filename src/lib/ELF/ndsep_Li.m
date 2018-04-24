function diimfp = ndsep_Li(osc,E0,depth,alpha)
%%
%{
   Calculates the normalised NDSEP (in eV^-1*A^-1) 
   for a given energy, angle and depth
   from solid to vacuum
   according to the Li algorithm Eq.(9)
   Y.C. Li et al. / Surface Science 589 (2005) 67-76.
%}
%%

x_in = zeros(length(osc.eloss),1);
w = osc.eloss;

for i=1:length(w)-1 
    qmin = sqrt(2*E0/h2ev)-sqrt(2*(E0/h2ev-w(i)/h2ev));
    qmax = sqrt(2*E0/h2ev)+sqrt(2*(E0/h2ev-w(i)/h2ev));   
    if (depth >= 0)
        sum4 = mu_s_v(osc,qmin/a0,qmax/a0,17,w(i),depth,alpha,E0,4);
        x_in(i) = 4*cosd(alpha)*sum4/pi^3;
    else
        sum1 = mu_s_v(osc,qmin/a0,qmax/a0,17,w(i),depth,alpha,E0,1); 
        sum2 = mu_s_v(osc,qmin/a0,qmax/a0,17,w(i),depth,alpha,E0,2);
        sum3 = mu_s_v(osc,qmin/a0,qmax/a0,17,w(i),depth,alpha,E0,3);
        x_in(i) = sum1/(pi*E0) - 2*cosd(alpha)*sum2/pi^3 + 4*cosd(alpha)*sum3/pi^3;
    end
end
diimfp = x_in ./ trapz(w,x_in);
plot(osc.eloss,diimfp)
hold on
plot(osc.eloss,x_in)
end