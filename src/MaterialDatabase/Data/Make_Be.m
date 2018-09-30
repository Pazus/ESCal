function Be = Make_Be
%% Be
Be.Mat = 'Be';
Be.M = 9.01218;
Be.Z = 4;
Be.Density = 1.848; %g/cm^3
Be.Density = Be.Density*10^-21/Be.M*6.022*10^23;
Be.NvTPP = 2;
Be.NvSGS = 2;
Be.Eg = 0;
Be.Ep = 18.44; %free-electron plasmon energy
Be.Ef = 14.3; %Fermi energy

E0 = [100:20:500 600:100:2500 2750:250:5000 5500:500:40000]';
Be.Elastic.x = zeros(numel(E0),1);
Be.Elastic.l_el = zeros(numel(E0),1);
Be.Elastic.l_tr = zeros(numel(E0),1);
Be.Elastic.x = E0;
Be.DECS.E0 = E0;

tic;
[data] = ElsepaRunner.RunElsepa(Be.Z,E0);
toc
Be.DECS.x = data(1).x;
for i = 1:numel(E0)
    Be.Elastic.l_el(i) = 1/data(i).sigma_el/Be.Density;
    Be.Elastic.l_tr(i) = 1/data(i).sigma_tr1/Be.Density;
    Be.DECS.y(:,i) = data(i).y/trapz(data(i).x,data(i).y);
end

Be.DIIMFP.E0 = [100;200;500;1000;2000;5000;10000;];
for i=1:numel(Be.DIIMFP.E0)
    WernerData = load([cd '/W_in/' Be.Mat num2str(Be.DIIMFP.E0(i)) '.diimfp']);
    Be.DIIMFP.y(:,i) = WernerData(:,3);
    if i==1
        Be.DIIMFP.x = WernerData(:,1);
    end
end

Be.XPS.E = [100; 200; 500; 1000; 1500; 2000; 3000; 4000; 5000];
Be.XPS.Shells = {'1S1/2','2S1/2'};
Be.XPS.EB = [115;9];

X = [.4680e3 .1773e3 .2972e2 .5539e1 .1870e1 .8315e0 .2525e0 .1053e0 .5267e-1;
    1.999 1.998 1.995 1.991 1.985 1.980 1.970 1.959 1.949;
    .160e0 .262e0 .463e0 .684e0 .855e0 .100e1 .124e1 .144e1 .162e1;
    -.405e-6 -.264e-6 .702e-7 .351e-6 .485e-6 .594e-6 .848e-6 .117e-5 .160e-5;
    .8391e2 .1884e2 .2054e1 .3027e0 .9275e-1 .3915e-1 .1127e-1 .4559e-2 .2239e-2;
    2.000 1.999 1.996 1.991 1.986 1.981 1.970 1.960 1.949;
    .112e0 .213e0 .427e0 .672e0 .846e0 .985e0 .123e1 .144e1 .162e1;
    -.782e-7 -.760e-7 -.693e-8 .225e-7 .687e-7 .159e-6 .346e-6 .620e-6 .961e-6;];

Be.XPS.PCS = X([0 4]+1, :)'*10^-7; %kbarn
Be.XPS.betta = X([0 4]+2, :)';
Be.XPS.gamma = X([0 4]+3, :)';
Be.XPS.delta = X([0 4]+4, :)';

