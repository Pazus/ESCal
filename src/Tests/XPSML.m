clear all

%% Material settings
Mat1 = XPSMaterial('Au',1); %Al anode as a default, you can set any shell here, it will be changed according to the shell settings  
Mat2 = XPSMaterial('Si',1);  
Mat1.CalculateDIIMFP; %DIIMFP from Werner data
Mat2.CalculateDIIMFP;

%% Geometry settings
theta0 = 0;
theta  = 45; phi=0;

%% Layers settings
Layers = [Layer(Mat1,1) Layer(Mat2)]; %two layers Au-Si
% Layers = [Layer(Mat1)]; %for one layer of Au
% Layers = [Layer(Mat1,3) Layer(Mat2,1) Layer(Mat1)]; %for three layers Au-Si-Au

%% Methods settings
% Methods = {'SLA'};
% Methods = {'SLA','NS'};
% Methods = {'NS','NS'};
% If you want to use NS you can skip this step

%% Shells settings
MS.Au={'4S1/2'};
MS.Si={'2S1/2'};
% MS.Au={'4S1/2','4P1/2','4P3/2'}; %for several shells
% MS.Si={'2S1/2','2P1/2','2P3/2'};

%% Run calculation
Rml = PESMultiLayer(Layers,MS);
%Rml = PESMultiLayer(Layers,MS,Methods);

Rml.theta0 = theta0; Rml.theta = theta; Rml.phi = phi;
Rml.N_in = 15;
Rml.energy_mesh_full = Mat1.Anode.PhotonEnergy-Mat1.BindingEnergy-50:0.5:Mat2.Anode.PhotonEnergy-Mat2.BindingEnergy +10;
Rml.sigma_gauss = 0.8;
Rml.Calculate;
% figure;
Rml.CalculateEnergyDistribution(theta,phi);
Rml.plotEnergyDistribution;

xlabel('Kinetic energy, eV');
ylabel('Intensity. Rel.un.');
