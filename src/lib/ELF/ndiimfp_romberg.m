function diimfp = ndiimfp_romberg(osc,E0)
%%
%{
   Calculates the normalised NDIIMFP (in eV^-1) 
   from data for the dielectric loss function, i.e. 
   the imaginary part of the reciprocal of the dielectric function.
   The Romberg algoritm is used.
%}
%%

x_in = zeros(length(osc.eloss),1);
% me = 9.10938356e-31;
% hbar = 6.582119514e-16;
w = osc.eloss;

for i=1:length(w) %-1
    qmin = sqrt(2*E0/h2ev)-sqrt(2*(E0/h2ev-w(i)/h2ev));
    qmax = sqrt(2*E0/h2ev)+sqrt(2*(E0/h2ev-w(i)/h2ev));
    sum = rombint(osc,qmin/a0,qmax/a0,15,w(i)); %Romberg intergration
%     osc.qtran = qmin/a0:0.01:qmax/a0;
%     osc.eloss = w(i);
%     ELF = eps_sum(osc);
%     sum = trapz(osc.qtran,ELF);
    x_in(i) = sum/(pi*E0);
    if ~mod(i,1000)
        YY = [num2str(i) ,' from ', num2str(length(w)-1)];
        disp(YY);
    end
end
diimfp = x_in ./ trapz(w,x_in);
%plot(osc.eloss,diimfp)

end