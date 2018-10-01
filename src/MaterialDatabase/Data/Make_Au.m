function Au = Make_Au()
%% Au
load('Au.mat');
Au.Mat = 'Au';
Au.M = 196.966569;
Au.Z = 79; %g/mol
Au.Density = 19.3; %g/cm^-3
Au.Density = Au.Density/Au.M*6.022*10^23*10^-21;
Au.NvTPP = 11;
Au.NvSGS = 11;
Au.Eg = 0;
Au.Ep = 29.92;
Au.Ef = 9.0;

Au.DIIMFP.E0 = [100;200;500;1000;2000;5000;10000;];
for i=1:numel(Au.DIIMFP.E0)
    WernerData = load([cd '\W_in\' Au.Mat num2str(Au.DIIMFP.E0(i)) '.diimfp']);
    Au.DIIMFP.y(:,i) = WernerData(:,3);
    if i==1
         Au.DIIMFP.x = WernerData(:,1);
    end
end

E0 = [100:20:500 600:100:2500 2750:250:5000 5500:500:40000]';
Au.Elastic.x = zeros(numel(E0),1);
Au.Elastic.l_el = zeros(numel(E0),1);
Au.Elastic.l_tr = zeros(numel(E0),1);
Au.Elastic.x = E0;
Au.DECS.E0 = E0;

tic;
[data] = ElsepaRunner.RunElsepa(Au.Z,E0);
toc
Au.DECS.x = data(1).x;
for i = 1:numel(E0)
    Au.Elastic.l_el(i) = 1/data(i).sigma_el/Au.Density;
    Au.Elastic.l_tr(i) = 1/data(i).sigma_tr1/Au.Density;
    Au.DECS.y(:,i) = data(i).y/trapz(data(i).x,data(i).y);
end

Au.XPS.Shells = {'4S1/2';'4P1/2';'4P3/2';'4D3/2';'4D5/2';'4F5/2';'4F7/2';'5S1/2';'5P1/2';'5P3/2';'5D3/2';'5D5/2';'6S1/2';};
Au.XPS.EB = [763;646;549;356;339;91;87;113;76;60;10;10;9;];
Au.XPS.E = [100; 200; 500; 1000; 1500; 2000; 3000; 4000; 5000;];
for i = 1:numel(Au.XPS.EB)
    Au.XPS.PCS(:,i) = S(:,i+3*(i-1))*10^-7;
    Au.XPS.betta(:,i) = S(:,i+3*(i-1)+1);
    Au.XPS.gamma(:,i) = S(:,i+3*(i-1)+2);
    Au.XPS.delta(:,i) = S(:,i+3*(i-1)+3);
end
