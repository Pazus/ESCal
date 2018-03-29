function optfit

l = load('hopgpal'); %load a file with Palik's data
hopg = l.cgraphite;

x=hopg(1:301,1);     %first column with energies
y_exp=hopg(1:301,4); %fourth column with Im[-1/epsilon]

Mat = XPSMaterial('C',1); %declare material properties

%% Oscillator parameters
osc.qtran = 0.09;     %momentum transfer should be zero, as we fit to optical data (except for the Mermin like models where it is almost zero just due to the calculation stuff)
osc.eloss = x;        %energy scale should be the same as experimental one
osc.model = 'Mermin'; %model:Drude, Drude-Lindhard, Mermin, Mermin_LL
osc.Ef = Mat.Ef;      %the Fermi energy

% The number of oscillators and their boundaries are the most important part of
% a successful fitting. Read Marteen Vos papers when setting all the
% parameters, e.g.
% Maarten Vos, Pedro L. Grande,
% How the choice of model dielectric function affects the calculated observables,
% Nuclear Instruments and Methods in Physics Research Section B: Beam Interactions with Materials and Atoms,
% V. 407, 2017, Pp. 97-109, https://doi.org/10.1016/j.nimb.2017.05.064.
% (http://www.sciencedirect.com/science/article/pii/S0168583X17306638)

osc.A = [0.1157 0.2004 0.1 0.1528 0.3480];   %amplitude
osc.G = [1.1258 8.3847 5   6.5495 8.3528];   %width (or damping)
osc.Om =[6.6639 5.3057 10  19.4540 25.4404]; %omega_pl ("position")

%lower boundaries
osc_min.A = zeros(size(osc.A));
osc_min.G = zeros(size(osc.A));
osc_min.Om =zeros(size(osc.A));
%upper boundaries
osc_max.A = Inf(size(osc.A));
osc_max.G = Inf(size(osc.A));
osc_max.Om =Inf(size(osc.A));

lb = structToVec(osc_min);
ub = structToVec(osc_max);

% alpha doesn't affect the fitting at all as it is used later when calculating
% ELF in order to add a dispersion, i.e. for q > 0
osc.alpha = 1;   %alpha is a constant between 1 (metals) and 0 (insulators)

%not necessarily to set, depends on the model
osc.u = 12;      %value responsible for the energy gap E_g (insulators and semiconductors), which is used only in Mermin_LL
osc.beps = 1;    %background dielectric constant. Is used only in Drude
osc.egap = 0;    %energy gap, is used only in Drude

%ionisation
osc.numion = 0;  %number of inner shells to be calculated
osc.ion = false; %calculate ionisation or not

%% Fit settings

function v = structToVec(s)
    v = [s.A, s.G, s.Om];
end

function s = vecToStruct(v)   
    s = osc;
    nA = length(osc.A);
    s.A = v(1:nA);
    s.G = v((1:nA)+nA);
    s.Om = v((1:nA)+nA*2);
end

function y = fit_func(x,xdata)
    o = vecToStruct(x);
    o.eloss = xdata;
    elf0 = eps_sum(o);
    elf = elf0(:,1);
    elf(1)=eps;
    y = elf;
    %plot(xdata,y);
end

%local minimisation
% {
options = optimoptions('lsqcurvefit','PlotFcn',@optimplotx,'Display','iter-detailed','Algorithm','levenberg-marquardt','UseParallel',true);
options.FunctionTolerance = 1.000000e-08;
[x_res,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@fit_func, structToVec(osc),x,y_exp,lb,ub,options);
%}

%global minimisation
%{
problem = createOptimProblem('lsqcurvefit','x0',structToVec(osc),'objective',@fit_func,...
    'lb',lb,'ub',ub,'xdata',x,'ydata',y_exp);
ms = MultiStart('PlotFcns',@gsplotbestf,'Display','iter','UseParallel',true,'XTolerance',1e-8);
[x_res,errormulti] = run(ms,problem,50)
%}

%% save and plot
figure;
an=vecToStruct(x_res);
save matnameosc.mat an

plot(x,y_exp,x,fit_func(x_res,x)); 

end