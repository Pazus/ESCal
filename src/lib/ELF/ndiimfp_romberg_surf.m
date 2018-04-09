function diimfp = ndiimfp_romberg_surf(osc,E0,theta)
%%
%{
   Calculates the normalised NDSEP (in eV^-1) 
   from data for the dielectric loss function, i.e. 
   the imaginary part of the reciprocal of the dielectric function.
   The Romberg algoritm is used.
%}
%%

x_in = zeros(length(osc.eloss),1);
w = osc.eloss;

for i=1:length(w) 
    qmin = sqrt(2*E0/h2ev)-sqrt(2*(E0/h2ev-w(i)/h2ev));
    qmax = sqrt(2*E0/h2ev)+sqrt(2*(E0/h2ev-w(i)/h2ev));
    sum = rombint_surf(osc,qmin/a0,qmax/a0,15,w(i),theta,E0); %Romberg intergration
    x_in(i) = sum/(pi*E0*cosd(theta));
    if ~mod(i,1000)
        YY = [num2str(i) ,' from ', num2str(length(w)-1)];
        disp(YY);
    end
end
diimfp = x_in ./ trapz(w,x_in);
%plot(osc.eloss,diimfp)
end