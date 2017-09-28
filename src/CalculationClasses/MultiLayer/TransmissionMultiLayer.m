classdef TransmissionMultiLayer<BaseMultiLayer
    
    properties
        %          z; %заменить все tau на z, вектор длиной в количество слоев
        theta                  = 0;
        phi                    = [];
        %         delta_max              = [];
        Fm;
        
    end
    
    methods
        function obj = TransmissionMultiLayer (Layers, CalculationMethod)
            
            obj = obj@BaseMultiLayer(Layers, CalculationMethod);
            
            obj.dE = abs(obj.Layers(1).Material.DIIMFP_E(2) - obj.Layers(1).Material.DIIMFP_E(1));
            obj.energy_mesh_full=obj.Layers(1).Material.DIIMFP_E;
            % Проверка на одинаковость E0  для всех слоёв
            %             for i_layer = 1:N_Layer
            %
            %             end
            %             for i = 1:obj.N_layer
            %                     obj.z(i) = Layers(i).thickness;
            %                     obj.Mat(i)=Layers(i).MaterialName;
            %             end
        end
        
        function Calculate(obj)
            for i_layer = 1:obj.N_Layer
                switch obj.CalculationMethod{i_layer}
                    case 'NS'
                        obj.ObjectsOfLayers{i_layer} = NSTransmition(obj.Layers(i_layer));
                    case 'SLA'
                        obj.ObjectsOfLayers{i_layer} = SLATransmition(obj.Layers(i_layer));
                    case 'SA'
                        obj.ObjectsOfLayers{i_layer} = SATransmition(obj.Layers(i_layer));
                    otherwise
                        error('Введите корректное название метода: NS, SLA или SA');
                end
            end
            obj.mu_mesh = obj.ObjectsOfLayers{1}.mu_mesh;
            
            obj.CalculateEachLayer;
            obj.Setup_energy_mesh_full;
            obj.RecalculateToEnergyDistribution;
            T = obj.FullEnergyDistribution{obj.N_Layer}.T(:,:,:,:)+obj.FullEnergyDistribution{obj.N_Layer}.L(:,:,:,:);
            w = sparse(1:obj.ObjectsOfLayers{1}.N, 1:obj.ObjectsOfLayers{1}.N, obj.ObjectsOfLayers{1}.mu_mesh_weights./obj.mu_mesh);
            M = size(T,4)-1;
            N_E =size(T,3);
            wT = zeros(size(obj.FullEnergyDistribution{1}.T));
            Tw = wT;
            R = wT;
            K_E =[];
            wR_cur = wT;
            Tw_cur = wT;
            for i_layer = (obj.N_Layer-1):-1:1
                if i_layer==obj.N_Layer-1
                    R = obj.FullEnergyDistribution{obj.N_Layer}.R(:,:,:,:);
                else
                    Rlocal = obj.FullEnergyDistribution{obj.i_layer+1}.R(:,:,:,:);
%                     R = obj.FullEnergyDistribution{obj.N_Layer}.R(:,:,:,:);
                    for m=0:M
                        for i=1:N_E
                            wT(:,:,i,m+1) = w*obj.FullEnergyDistribution{i_layer+1}.T(:,:,i,m+1) + obj.FullEnergyDistribution{i_layer+1}.L(:,:,i,m+1);
                            Tw(:,:,i,m+1) = obj.FullEnergyDistribution{i_layer+1}.T(:,:,i,m+1)*w + obj.FullEnergyDistribution{i_layer+1}.L(:,:,i,m+1);
                        end
                    end
                    R = obj.DoubleLayer(Rlocal,Tw,R,wT,K_E);
                end
                
                for m=0:M
                    for i=1:N_E
                        wR_cur(:,:,i,m+1) = obj.FullEnergyDistribution{i_layer}.R(:,:,i,m+1)*w;
                        Tw_cur(:,:,i,m+1) = obj.FullEnergyDistribution{i_layer}.T(:,:,i,m+1)*w + obj.FullEnergyDistribution{i_layer}.L(:,:,i,m+1);
                        T     (:,:,i,m+1) = w*T(:,:,i,m+1);
                    end
                end
                K_E = length(obj.energy_mesh_full)-find(obj.energy_mesh_full==obj.ObjectsOfLayers{1}.Material.E0); % !!!исправить find!!! Число элементов после E0
                C = repmat(eye(size(obj.FullEnergyDistribution{1}.T,1)),1,1,N_E,M+1) + conv3d(wR_cur,R,'right',K_E)*obj.dE;
                T = obj.DoubleLayer(0,Tw_cur,C,T,K_E);
            end
            obj.Fm = T;
            obj.IsCalculated = true;
        end
        
        function CalculateEnergyDistribution(obj,theta,phi)
            if nargin < 3; phi = 0; end; %??????
            obj.CalculateEnergyDistribution_E0(theta,phi);
            convGauss (obj,obj.sigma_gauss);
            convLDS (obj,obj.sigma_LDS, obj.alpha_LDS);
        end
        
    end % methods
end % classdef

