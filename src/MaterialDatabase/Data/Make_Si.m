function Si = Make_Si
%% Si
load('Si.mat');
Si.Mat = 'Si';
Si.M = 28.0855;
Si.Z = 14;
Si.Density = 2.33; %g/cm^3
Si.Density = Si.Density*10^-21/Si.M*6.022*10^23;
Si.NvTPP = 4;
Si.NvSGS = 4;
Si.Eg = 1.1;
Si.Ep = 16.59;
Si.Ef = 12.5;

Si.DIIMFP.E0 = [100;200;500;1000;2000;5000;10000;];
for i=1:numel(Si.DIIMFP.E0)
    WernerData = load([cd '/W_in/' Si.Mat num2str(Si.DIIMFP.E0(i)) '.diimfp']);
    Si.DIIMFP.y(:,i) = WernerData(:,3);
    if i==1
         Si.DIIMFP.x = WernerData(:,1);
    end
end

E0 = [100:20:500 600:100:2500 2750:250:5000 5500:500:40000]';
Si.Elastic.x = zeros(numel(E0),1);
Si.Elastic.l_el = zeros(numel(E0),1);
Si.Elastic.l_tr = zeros(numel(E0),1);
Si.Elastic.x = E0;
Si.DECS.E0 = E0;

tic;
[data] = ElsepaRunner.RunElsepa(Si.Z,E0);
toc
Si.DECS.x = data(1).x;
for i = 1:numel(E0)
    Si.Elastic.l_el(i) = 1/data(i).sigma_el/Si.Density;
    Si.Elastic.l_tr(i) = 1/data(i).sigma_tr1/Si.Density;
    Si.DECS.y(:,i) = data(i).y/trapz(data(i).x,data(i).y);
end

Si.XPS.Shells = {'2S1/2';'2P1/2';'2P3/2';'3S1/2';'3P1/2';};
Si.XPS.EB = [154;104;104;13;8;];
Si.XPS.E = [100; 200; 500; 1000; 1500; 2000; 3000; 4000; 5000;];
S(:,4)=S(:,4)*10^-5;
S(:,16)=S(:,16)*10^-5;
for i = 1:numel(Si.XPS.EB)
    Si.XPS.PCS(:,i) = S(:,i+3*(i-1))*10^-7;
    Si.XPS.betta(:,i) = S(:,i+3*(i-1)+1);
    Si.XPS.gamma(:,i) = S(:,i+3*(i-1)+2);
    Si.XPS.delta(:,i) = S(:,i+3*(i-1)+3);
end