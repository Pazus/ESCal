function Mo = Make_Mo()
%% Cu
load('Mo.mat');
Mo.Mat = 'Mo';
Mo.M = 95.94;
Mo.Z = 42; %g/mol
Mo.Density = 10.28; %g/cm^-3
Mo.Density = Mo.Density/Mo.M*6.022*10^23*10^-21;
Mo.NvTPP = 6;
Mo.NvSGS = 6;
Mo.Eg = 0;
Mo.Ep = 23.09;
Mo.Ef = 6.5;

Mo.DIIMFP.E0 = [200;500;1000;2000;5000;10000;];
for i=1:numel(Mo.DIIMFP.E0)
    WernerData = load([cd '\W_in\' Mo.Mat num2str(Mo.DIIMFP.E0(i)) '.diimfp']);
    Mo.DIIMFP.y(:,i) = WernerData(:,3);
    if i==1
         Mo.DIIMFP.x = WernerData(:,1);
    end
end

E0 = [100:20:500 600:100:2500 2750:250:5000 5500:500:40000]';
Mo.Elastic.x = zeros(numel(E0),1);
Mo.Elastic.l_el = zeros(numel(E0),1);
Mo.Elastic.l_tr = zeros(numel(E0),1);
Mo.Elastic.x = E0;
Mo.DECS.E0 = E0;

tic;
[data] = ElsepaRunner.RunElsepa(Mo.Z,E0);
toc
Mo.DECS.x = data(1).x;
for i = 1:numel(E0)
    Mo.Elastic.l_el(i) = 1/data(i).sigma_el/Mo.Density;
    Mo.Elastic.l_tr(i) = 1/data(i).sigma_tr1/Mo.Density;
    Mo.DECS.y(:,i) = data(i).y/trapz(data(i).x,data(i).y);
end

Mo.XPS.Shells = {'1S1/2';'2S1/2';'2P1/2';'2P3/2';'3S1/2';'3P1/2';'3P3/2';'3D3/2';'3D5/2';'4S1/2';'4P1/2';'4P3/2';'4D3/2';'4D5/2';'5S1/2';};
Mo.XPS.EB = [20006;2872;2631;2526;511;416;399;236;234;68;44;41;8;8;7;];
Mo.XPS.E = [100; 200; 500; 1000; 1500; 2000; 3000; 4000; 5000;];
% S(:,4)=S(:,4)*10^-4;
% S(:,16)=S(:,16)*10^-4;
% S(:,36)=S(:,36)*10^-4;
for i = 1:numel(Mo.XPS.EB)
    Mo.XPS.PCS(:,i) = S(1+(i-1)*4,:)'*10^-7;
    Mo.XPS.betta(:,i) = S(2+(i-1)*4,:);
    Mo.XPS.gamma(:,i) = S(3+(i-1)*4,:);
    Mo.XPS.delta(:,i) = S(4+(i-1)*4,:);
end