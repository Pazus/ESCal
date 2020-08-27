function Li = Make_Li
%% Be
Li.Mat = 'Li';
Li.M = 6.94; % [6.938, 6.997] 
Li.Z = 3;
Li.Density = 0.534; %g/cm^3
Li.Density = Li.Density*10^-21/Li.M*6.022*10^23;
Li.NvTPP = 1;
Li.NvSGS = 1;
Li.Eg = 0;
Li.Ep = 7.99; %free-electron plasmon energy
Li.Ef = 4.74; %Fermi energy

E0 = [100:20:500 600:100:2500 2750:250:5000 5500:500:40000]';
Li.Elastic.x = zeros(numel(E0),1);
Li.Elastic.l_el = zeros(numel(E0),1);
Li.Elastic.l_tr = zeros(numel(E0),1);
Li.Elastic.x = E0;
Li.DECS.E0 = E0;

tic;
[data] = ElsepaRunner.RunElsepa(Li.Z,E0);
toc
Li.DECS.x = data(1).x;
for i = 1:numel(E0)
    Li.Elastic.l_el(i) = 1/data(i).sigma_el/Li.Density;
    Li.Elastic.l_tr(i) = 1/data(i).sigma_tr1/Li.Density;
    Li.DECS.y(:,i) = data(i).y/trapz(data(i).x,data(i).y);
end

% Li.DIIMFP.E0 = [100;200;500;1000;2000;5000;10000;];
% for i=1:numel(Li.DIIMFP.E0)
%     WernerData = load([cd '\W_in\' Li.Mat num2str(Li.DIIMFP.E0(i)) '.diimfp']);
%     Li.DIIMFP.y(:,i) = WernerData(:,3);
%     if i==1
%         Li.DIIMFP.x = WernerData(:,1);
%     end
% end

Li.XPS.E = [100; 200; 500; 1000; 1500; 2000; 3000; 4000; 5000];
Li.XPS.Shells = {'1S1/2','2S1/2'};
Li.XPS.EB = [58;5];

X = [.3658e3 .1044e3 .1276e2 .1996e1 .6168e0 .2601e0 .7406e-1 .2976e-1 .1451e-1;
    1.999 1.998 1.996 1.990 1.985 1.980 1.969 1.959 1.948;
    .177e0 .276e0 .471e0 .699e0 .872e0 .102e1 .126e1 .146e1 .163e1;
    -.121e-6 .180e-8 .195e-6 .262e-6 .302e-6 .328e-6 .494e-6 .700e-6 .106e-5;
    .2905e2 .6171e1 .5399e0 .07194 .02094 .008490 .002347 .0009237 .0004462;
    2.000 1.999 1.996 1.991 1.985 1.980 1.969 1.959 1.947;
    .141e0 .238e0 .468e0 .683e0 .879e0 .101e1 .127e1 .145e1 .164e1;
    -0.534e-7 .228e-7 -.469e-8 .166e-7 .524e-7 .127e-6 .201e-6 .603e-6 .583e-6;];

Li.XPS.PCS = X([0 4]+1, :)'*10^-7; %kbarn
Li.XPS.betta = X([0 4]+2, :)';
Li.XPS.gamma = X([0 4]+3, :)';
Li.XPS.delta = X([0 4]+4, :)';

