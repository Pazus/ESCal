clear;

osc.A = [64 6 6.5 5.5 4 55 42 172 80 240 90 85 200 500 664];
osc.G = [0.03 0.3 0.65 0.7 0.7 2.6 4.76 10.18 8 32 30 30 25 65 160];
osc.Om = [0 0.3 2.5 3.1 3.7 5.05 8.93 14.74 25.6 40 55 65 83 120 200];

osc.model = 'Drude';
osc.alpha = 1; 
osc.beps = 1.05;
osc.Ef = 7;
osc.egap = 0;
osc.numion = 0;
osc.ion = false;

E0 = 500;
depth = 0;
alpha = 0;
an.eloss = 0:.5:E0;

[diimfp,dsep] = ndiimfp_Li(an,E0,depth,alpha,12);

figure;
plot(an.eloss,diimfp)
hold on
plot(an.eloss,dsep)
legend('ndiimfp','ndsep')
xlim([0 100])