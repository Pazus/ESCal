function Al = Make_Al()
%% Al
load('Al.mat');
Al.Mat = 'Al';
Al.M = 26.9815386;
Al.Z = 13;
Al.Density = 2.70; %g/cm^3
Al.Density = Al.Density*10^-21/Al.M*6.022*10^23;
Al.NvTPP = 3;
Al.NvSGS = 3;
Al.Eg = 0;
Al.Ep = 15.78;
Al.Ef = 11.2;

Al.DIIMFP.E0 = [100;200;500;1000;2000;5000;10000;];
for i=1:numel(Al.DIIMFP.E0)
    WernerData = load([cd '/W_in/' Al.Mat num2str(Al.DIIMFP.E0(i)) '.diimfp']);
    Al.DIIMFP.y(:,i) = WernerData(:,3);
    if i==1
         Al.DIIMFP.x = WernerData(:,1);
    end
end

E0 = [100:20:500 600:100:2500 2750:250:5000 5500:500:40000]';
Al.Elastic.x = zeros(numel(E0),1);
Al.Elastic.l_el = zeros(numel(E0),1);
Al.Elastic.l_tr = zeros(numel(E0),1);
Al.Elastic.x = E0;
Al.DECS.E0 = E0;

tic;
[data] = ElsepaRunner.RunElsepa(Al.Z,E0);
toc
Al.DECS.x = data(1).x;
for i = 1:numel(E0)
    Al.Elastic.l_el(i) = 1/data(i).sigma_el/Al.Density;
    Al.Elastic.l_tr(i) = 1/data(i).sigma_tr1/Al.Density;
    Al.DECS.y(:,i) = data(i).y/trapz(data(i).x,data(i).y);
end

Al.XPS.Shells = {'2S1/2';'2P1/2';'2P3/2';'3S1/2';'3P1/2';};
Al.XPS.EB = [121;77;77;10;6;];
Al.XPS.E = [100; 200; 500; 1000; 1500; 2000; 3000; 4000; 5000;];
S(:,4)=S(:,4)*10^-5;
S(:,16)=S(:,16)*10^-5;
for i = 1:numel(Al.XPS.EB)
    Al.XPS.PCS(:,i) = S(:,i+3*(i-1))*10^-7;
    Al.XPS.betta(:,i) = S(:,i+3*(i-1)+1);
    Al.XPS.gamma(:,i) = S(:,i+3*(i-1)+2);
    Al.XPS.delta(:,i) = S(:,i+3*(i-1)+3);
end