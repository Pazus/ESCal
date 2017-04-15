classdef XPSMultiLayer < handle
% x_in.delta - сетка по потерям энергии для индикатрисы
% x_in.I - матрица индикатрис
% x_in.I(i,:) - сечение для i-го слоя, начиная с полубесконечности (i=1)

% TODO
% На вход передавать массив Layer
% Каждый Layer должен быть из одного и того же материала
% Данное ограничени необходимо контролировать проверкой в конструкторе
% (это временное решение, потом ограничени будет снято)
    properties
        Mat@XPSMaterial;
        Q
        mesh_E
        x_in
        z
        N_layer                
    %end
      
    %properties (SetObservable = true, AbortSet = true)
        N_in                   = 30;
        theta0                 = 0;
        theta                  = 0;
        phi                    = 0;
        sigma_gauss            = [];
        sigma_DS               = [];
        alpha_DS               = []
        dE                     = 0.5;
        deltaBE_exp            = [];
        nShellsCalc            = [];
    end 
    
    properties (SetObservable = true, AbortSet = true)
        mesh_E_global = [];
    end
    
    
    methods
        
        function obj = XPSMultiLayer(Mat,x_in,z)

            if nargin<3;
                obj.N_layer = 1;
            else
                obj.N_layer = 1 + length(z);
                obj.z=z;                
            end
                        
            obj.Mat=Mat;
            obj.x_in=x_in;

            addlistener(obj,'mesh_E_global','PostSet', @obj.Recalculate_dE);
        end;
        
        function Recalculate_dE(obj,varargin)
           obj.dE = abs(obj.mesh_E_global(2)-obj.mesh_E_global(1)); 
        end        

        function Calculate(obj)
            
            function y = calcDIIMFP(x0,I_in,mesh_E,mesh_E_global)
                I_in_temp = interp1(x0,I_in,mesh_E);
                I_in_temp(isnan(I_in_temp))=0;

                y = fliplr(I_in_temp/trapz(mesh_E,I_in_temp));
                y = interp1(mesh_E,y,mesh_E_global);

                y(isnan(y))=0;
            end
                N_layer = obj.N_layer;
                theta  = obj.theta;
                phi    = obj.phi;
                theta0 = obj.theta0;
                N_in   = obj.N_in;
                dE     = obj.dE;
                Mat    = obj.Mat;
                mesh_E_global = obj.mesh_E_global;
                nShellsCalc = obj.nShellsCalc;
                deltaBE_exp = obj.deltaBE_exp;
                sigma_gauss = obj.sigma_gauss;
                sigma_DS = obj.sigma_DS;
                alpha_DS = obj.alpha_DS;
                x_in = obj.x_in;
                if obj.N_layer>1 z = obj.z; end
                if isempty(nShellsCalc) nShellsCalc = 1:Mat.ShellsCount;end
                if isempty(sigma_gauss) sigma_gauss = zeros(1,Mat.ShellsCount);end
                if isempty(sigma_DS) sigma_DS = zeros(1,Mat.ShellsCount);end
                if isempty(deltaBE_exp) deltaBE_exp = zeros(1,Mat.ShellsCount);end
                if isempty(mesh_E_global) mesh_E_global = 0:dE:Mat.Anode.PhotonEnergy;end
                Q_global = zeros(1,length(mesh_E_global));
       
                for n = nShellsCalc
                Mat.SetShell(n);
                Mat.BindingEnergy = Mat.BindingEnergy + deltaBE_exp(n);
                mesh_E_temp = 0:dE:Mat.E0;
                DIIMFP = zeros(N_layer,length(mesh_E_global));
                
                for i_layer = 1:N_layer
                    DIIMFP(i_layer,:) = calcDIIMFP(x_in.delta,x_in.I(i_layer,:),mesh_E_temp,mesh_E_global);
                end
                
%                 y_b = calcDIIMFP(x_in.delta,x_in.I(1,:),mesh_E_temp,mesh_E_global);
%                 y_s = calcDIIMFP(x_in.delta,x_in.I(2,:),mesh_E_temp,mesh_E_global);
                
                %% Q Bulk calculation 
                Mat.SetManualDIIMFP(mesh_E_global', DIIMFP(1,:)');
                Qc = NSXPS(Layer(Mat)); Qc.N_in = N_in; Qc.theta0 = theta0; 
                tic ;
                Qc.Calculate;
                toc 
                Qc.CalculateEnergyDistribution(theta,phi);
                Q_b = sum(Qc.EnergyDistribution');
    
                %% Surface
                
                for i_layer = 2:N_layer
                    Mat.SetManualDIIMFP(mesh_E_global', DIIMFP(i_layer,:)');
                    Q_s = SLAXPS(Layer(Mat,z(i_layer-1))); Q_s.N_in = N_in; Q_s.theta0 = theta0; 
                    Q_s.Calculate;
                    Q_s.CalculateEnergyDistribution(theta,phi)
                    Qs = sum(Q_s.EnergyDistribution');
                    Q_s.CalculationResultPropertyName = 'Tm';
                    Q_s.theta0 = theta; Q_s.Calculate;
                    Q_s.CalculateEnergyDistribution(theta,phi)
                    Ts = sum(Q_s.EnergyDistribution');
    
                    QbTs = conv_my (Q_b, Ts, dE);
                    Q_b = QbTs + Qs;
                end %по слоям
                
                %% Свёртка с гаусссом
                Q_res = zeros(size(Q_b));
                if sigma_gauss(n)==0
                    Q_res = Q_b;
                else
                    GaussF = normpdf(-3*sigma_gauss(n):dE:3*sigma_gauss(n), 0, sigma_gauss(n));
                    Q_res = conv_my(Q_b, GaussF,dE,'same');    
                end 
                
                if ~(sigma_DS(n)==0 || isempty(alpha_DS) || isempty(sigma_DS))
                    mesh_LorentzDS = -Mat.BindingEnergy:dE:Mat.BindingEnergy;
                    LorentzDS = gamma(1-alpha_DS(n))*cos(pi*alpha_DS(n)/2 + (1-alpha_DS(n))*atan(mesh_LorentzDS...
                               ./(sigma_DS(n)/2)))./(mesh_LorentzDS.^2 + (sigma_DS(n)/2).^2).^((1-alpha_DS(n))./2);
                    LorentzDS = LorentzDS/trapz(mesh_LorentzDS,LorentzDS);
                    Q_res = conv_my(Q_res, LorentzDS,dE,'same');
                end 
                Q_global = Q_global + Q_res;
                end % по оболочкам
            obj.Q = Q_global;
            obj.mesh_E = mesh_E_global;
        end
        
        function plotSpectrum(obj)
            plot(obj.mesh_E,obj.Q);
            xlabel('{\itE}, эВ')
            ylabel('Интенсивность  {\itI}, отн.ед.')
%             legend('Расчёт','Эксперимент' )
        end

    end
end

