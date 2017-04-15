function Ta = Make_Ta()
%% Cu
load('Ta.mat');
Ta.Mat = 'Ta';
Ta.M = 180.9479;
Ta.Z = 73; %g/mol
Ta.Density = 16.65; %g/cm^-3
Ta.Density = Ta.Density/Ta.M*6.022*10^23*10^-21;
Ta.NvTPP = 5;
Ta.NvSGS = 5;
Ta.Eg = 0;
Ta.Ep = 19.53;
Ta.Ef = 8.4;

Ta.DIIMFP.E0 = [100;200;500;1000;2000;5000;10000;];
for i=1:numel(Ta.DIIMFP.E0)
    WernerData = load([cd '\W_in\' Ta.Mat num2str(Ta.DIIMFP.E0(i)) '.diimfp']);
    Ta.DIIMFP.y(:,i) = WernerData(:,3);
    if i==1
         Ta.DIIMFP.x = WernerData(:,1);
    end
end

E0 = [100:20:500 600:100:2500 2750:250:5000 5500:500:40000]';
Ta.Elastic.x = zeros(numel(E0),1);
Ta.Elastic.l_el = zeros(numel(E0),1);
Ta.Elastic.l_tr = zeros(numel(E0),1);
Ta.Elastic.x = E0;
Ta.DECS.E0 = E0;

tic;
[data] = ElsepaRunner.RunElsepa(Ta.Z,E0);
toc
Ta.DECS.x = data(1).x;
for i = 1:numel(E0)
    Ta.Elastic.l_el(i) = 1/data(i).sigma_el/Ta.Density;
    Ta.Elastic.l_tr(i) = 1/data(i).sigma_tr1/Ta.Density;
    Ta.DECS.y(:,i) = data(i).y/trapz(data(i).x,data(i).y);
end

Ta.XPS.Shells = {'3D3/2';'3D5/2';'4S1/2';'4P1/2';'4P3/2';'4D3/2';'4D5/2';'4F5/2';'4F7/2';'5S1/2';'5P1/2';'5P3/2';'5D3/2';'6S1/2';};
Ta.XPS.EB = [1796;1737;570;469;407;245;232;30;28;74;47;38;8;8];
Ta.XPS.E = [100; 200; 500; 1000; 1500; 2000; 3000; 4000; 5000;];
% S(:,4)=S(:,4)*10^-4;
% S(:,16)=S(:,16)*10^-4;
% S(:,36)=S(:,36)*10^-4;
for i = 1:numel(Ta.XPS.EB)
    Ta.XPS.PCS(:,i) = S(1+(i-1)*4,:)'*10^-7;
    Ta.XPS.betta(:,i) = S(2+(i-1)*4,:);
    Ta.XPS.gamma(:,i) = S(3+(i-1)*4,:);
    Ta.XPS.delta(:,i) = S(4+(i-1)*4,:);
end