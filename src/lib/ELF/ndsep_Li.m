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
    if (depth > 0) % outside
        sum4 = mu_s_v(osc,qmin/a0,qmax/a0,15,w(i),depth,alpha,E0,4);
        x_in(i) = 4*cosd(alpha)*sum4/pi^3/a0/h2ev;
    elseif (depth < 0) % inside
        %sum1 = mu_s_v(osc,qmin/a0,qmax/a0,15,w(i),depth,alpha,E0,1); 
        %sum2 = mu_s_v(osc,qmin/a0,qmax/a0,15,w(i),depth,alpha,E0,2);
        sum3 = mu_s_v(osc,qmin/a0,qmax/a0,15,w(i),depth,alpha,E0,3);
        %clear_bulk(i) = sum1/(pi*E0)/a0/h2ev;
        %bulks_v(i) = (-2)*cosd(alpha)*sum2/pi^3/a0/h2ev;
        x_in(i) = 4*cosd(alpha)*sum3/pi^3/a0/h2ev;
    elseif (depth == 0)
        sum4 = mu_s_v(osc,qmin/a0,qmax/a0,15,w(i),depth,alpha,E0,4);
        sum3 = mu_s_v(osc,qmin/a0,qmax/a0,15,w(i),depth,alpha,E0,3);
        x_in(i) = 0.5*4*cosd(alpha)/pi^3/a0/h2ev*(sum4 + sum3);
    end
end
diimfp = x_in ./ trapz(w,x_in);
plot(w,diimfp)
hold on
plot(w,x_in)
end