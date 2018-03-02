clear;
osc = struct;

osc.A = [220 11];
osc.G = [01.4 25];
osc.Om = [4.5 26];

osc.alpha = 500;
osc.alpha = 1500;
osc.beps = 1;
osc.model = 'Drude';
osc.u = 12;
osc.Ef = 11.2;
osc.qtran = 0:10;
osc.eloss = 0:0.5:1000;
osc.egap = 0;

tic;
elf = eps_sum(osc);
toc
tic;
elf2 = eps_sum(osc);
toc

norm(elf-elf2) % ~1e-16
norm(elf)
norm(elf2) % >1

%surf(elf-elf2)
%semilogy(osc.eloss,elf)