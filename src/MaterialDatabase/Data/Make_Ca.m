function Ca = Make_Ca()
%% Ca
% load('Ña.mat');
Ca.Mat = 'Ña';
Ca.M = 40.078;
Ca.Z = 20; %g/mol
Ca.Density = 1.55; %g/cm^-3
Ca.Density = Ca.Density/Ca.M*6.022*10^23*10^-21;
% TPP2M doean't have cross sections
Ca.NvTPP = 2;
Ca.NvSGS = 2;
Ca.Eg = 0;
Ca.Ep = 8;
Ca.Ef = 4;

% Ca.DIIMFP.E0 = [100;200;500;1000;2000;5000;10000;];
% for i=1:numel(Ca.DIIMFP.E0)
%     WernerData = load([cd '\W_in\' Ca.Mat num2str(Ca.DIIMFP.E0(i)) '.diimfp']);
%     Ca.DIIMFP.y(:,i) = WernerData(:,3);
%     if i==1
%          Ca.DIIMFP.x = WernerData(:,1);
%     end
% end

E0 = [100:20:500 600:100:2500 2750:250:5000 5500:500:40000]';
Ca.Elastic.x = zeros(numel(E0),1);
Ca.Elastic.l_el = zeros(numel(E0),1);
Ca.Elastic.l_tr = zeros(numel(E0),1);
Ca.Elastic.x = E0;
Ca.DECS.E0 = E0;

tic;
[data] = ElsepaRunner.RunElsepa(Ca.Z,E0);
toc
Ca.DECS.x = data(1).x;
for i = 1:numel(E0)
    Ca.Elastic.l_el(i) = 1/data(i).sigma_el/Ca.Density;
    Ca.Elastic.l_tr(i) = 1/data(i).sigma_tr1/Ca.Density;
    Ca.DECS.y(:,i) = data(i).y/trapz(data(i).x,data(i).y);
end

% Ca.XPS.Shells = {'3D3/2';'3D5/2';'4S1/2';'4P1/2';'4P3/2';'4D3/2';'4D5/2';'4F5/2';'4F7/2';'5S1/2';'5P1/2';'5P3/2';'5D3/2';'6S1/2';};
% Ca.XPS.EB = [1796;1737;570;469;407;245;232;30;28;74;47;38;8;8];
% Ca.XPS.E = [100; 200; 500; 1000; 1500; 2000; 3000; 4000; 5000;];
% % S(:,4)=S(:,4)*10^-4;
% % S(:,16)=S(:,16)*10^-4;
% % S(:,36)=S(:,36)*10^-4;
% for i = 1:numel(Ca.XPS.EB)
%     Ca.XPS.PCS(:,i) = S(1+(i-1)*4,:)'*10^-7;
%     Ca.XPS.betta(:,i) = S(2+(i-1)*4,:);
%     Ca.XPS.gamma(:,i) = S(3+(i-1)*4,:);
%     Ca.XPS.delta(:,i) = S(4+(i-1)*4,:);
% end