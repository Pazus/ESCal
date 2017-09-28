clear all

%% Material settings
Mat1 = Material('Au',5000);  
Mat2 = Material('Si',5000);  
Mat1.CalculateDIIMFP; %DIIMFP from Werner data
Mat2.CalculateDIIMFP;

%% Geometry settings
theta0 = 0;
theta  = 45; phi=0;

%% Layers settings
Layers = [Layer(Mat1,3) Layer(Mat2)]; %two layers Au-Si
% Layers = [Layer(Mat1)]; %for one layer of Au
% Layers = [Layer(Mat1,3) Layer(Mat2,1) Layer(Mat1)]; %for three layers Au-Si-Au

%% Methods settings
% Methods = {'SLA'};
% Methods = {'SLA','NS'};
% Methods = {'NS','NS'};
% If you want to use NS you can skip this step

%% Run calculation
Rml = ReflectionMultiLayer(Layers);
%Rml = ReflectionMultiLayer(Layers,Methods);
% Rml = TransmissionMultiLayer(Layers,Methods);

Rml.theta0 = theta0; Rml.theta = theta; Rml.phi = phi;
Rml.N_in = 15;
Rml.energy_mesh_full = Mat1.E0-50:0.5:Mat1.E0+10;
Rml.sigma_gauss = 0.8;
Rml.Calculate;
% figure;
Rml.CalculateEnergyDistribution(theta,phi);
Rml.plotEnergyDistribution;

xlabel('Kinetic energy, eV');
ylabel('Intensity. Rel.un.');
