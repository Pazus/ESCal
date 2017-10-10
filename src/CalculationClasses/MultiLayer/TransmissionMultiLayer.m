classdef TransmissionMultiLayer_v2<BaseMultiLayer
    
    properties
        %          z; %заменить все tau на z, вектор длиной в количество слоев
        theta                  = 0;
        phi                    = [];
        %         delta_max              = [];
        Fm;
        
    end
    
    methods
        function obj = TransmissionMultiLayer_v2 (Layers, CalculationMethod)
            if nargin==1 || isempty(CalculationMethod);
                CalculationMethod = {'NS'};             
            end
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
            obj.CalculateEachLayer;
            obj.mu_mesh = obj.ObjectsOfLayers{1}.mu_mesh;
            
            obj.Setup_energy_mesh_full;
            obj.RecalculateToEnergyDistribution;
            w = sparse(1:obj.ObjectsOfLayers{1}.N, 1:obj.ObjectsOfLayers{1}.N, obj.ObjectsOfLayers{1}.mu_mesh_weights./obj.mu_mesh);
            wT = zeros(size(obj.FullEnergyDistribution{1}.T));
            Tw = wT;
            R = wT; K_E =[];
            Rw_cur = wT;
            Tw_cur = wT;
            T    = wT;
            Temp = obj.FullEnergyDistribution{obj.N_Layer}.T(:,:,:,:);
            M = size(Temp,4)-1;
            N_E =size(Temp,3); clear Temp;
            if strcmp(obj.CalculationMethod{obj.N_Layer},'SLA')
                T = obj.FullEnergyDistribution{obj.N_Layer}.T(:,:,:,:);
            else
                for m=0:M
                    for i=1:N_E
                        if obj.N_Layer==1
                            T (:,:,i,m+1)= obj.FullEnergyDistribution{obj.N_Layer}.T(:,:,i,m+1)+obj.FullEnergyDistribution{obj.N_Layer}.L(:,:,i,m+1);
                        else
                            T (:,:,i,m+1)= w*obj.FullEnergyDistribution{obj.N_Layer}.T(:,:,i,m+1)+obj.FullEnergyDistribution{obj.N_Layer}.L(:,:,i,m+1);
                        end
                    end
                end
            end
            
            for i_layer = (obj.N_Layer-1):-1:1
                %% расчёт вспомогательной функции отражения
                if i_layer==obj.N_Layer-1
                    R = obj.FullEnergyDistribution{obj.N_Layer}.R(:,:,:,:);
                else
                    Rlocal = obj.FullEnergyDistribution{i_layer+1}.R(:,:,:,:);
                    if strcmp(obj.CalculationMethod{i_layer+1},'SLA')
                        wT = obj.FullEnergyDistribution{i_layer+1}.T;
                        Tw = obj.FullEnergyDistribution{i_layer+1}.T;
                    else
                        for m=0:M
                            for i=1:N_E
                                wT(:,:,i,m+1) = w*obj.FullEnergyDistribution{i_layer+1}.T(:,:,i,m+1) + obj.FullEnergyDistribution{i_layer+1}.L(:,:,i,m+1);
                                Tw(:,:,i,m+1) = obj.FullEnergyDistribution{i_layer+1}.T(:,:,i,m+1)*w + obj.FullEnergyDistribution{i_layer+1}.L(:,:,i,m+1);
                            end
                        end
                    end
                    R = obj.DoubleLayer(Rlocal,Tw,R,wT,K_E);
                end
                %% расчёт пропускания
                if strcmp(obj.CalculationMethod{i_layer},'SLA')
                    Rw_cur = obj.FullEnergyDistribution{i_layer}.R;
                    Tw_cur = obj.FullEnergyDistribution{i_layer}.T;
                else
                    for m=0:M
                        for i=1:N_E
                            if i_layer==obj.N_Layer-1
                                Rw_cur(:,:,i,m+1) = obj.FullEnergyDistribution{i_layer}.R(:,:,i,m+1);
                            else
                                Rw_cur(:,:,i,m+1) = obj.FullEnergyDistribution{i_layer}.R(:,:,i,m+1)*w;
                            end
                            Tw_cur(:,:,i,m+1) = obj.FullEnergyDistribution{i_layer}.T(:,:,i,m+1)*w + obj.FullEnergyDistribution{i_layer}.L(:,:,i,m+1);
%                             if i_layer~=obj.N_Layer
                            R (:,:,i,m+1) = R(:,:,i,m+1)*w;
                                %    T (:,:,i,m+1) = w*T(:,:,i,m+1);
%                             end
                        end
                    end
                end
                K_E = length(obj.energy_mesh_full)-find(obj.energy_mesh_full==obj.ObjectsOfLayers{1}.Material.E0); % !!!исправить find!!! Число элементов после E0
                %{
                % вариант 1
                C = repmat(eye(size(obj.FullEnergyDistribution{1}.T,1)),1,1,N_E,M+1) + conv3d(R,Rw_cur,'right',K_E)*obj.dE;
                T = obj.DoubleLayer(0,Tw_cur,C,T,K_E);
                %}
                % {
                % вариант 2
                if i_layer==obj.N_Layer-1
                    % T_cur x T = (T1 + L1) x (T2 + L2) = T1 x (wT2 + L2) + L1 x (T2 + L2)
                    for m=0:M
                        for i=1:N_E
                            wT2L (:,:,i,m+1)= w*obj.FullEnergyDistribution{obj.N_Layer}.T(:,:,i,m+1)+obj.FullEnergyDistribution{obj.N_Layer}.L(:,:,i,m+1);
                            T2L (:,:,i,m+1)= obj.FullEnergyDistribution{obj.N_Layer}.T(:,:,i,m+1)+obj.FullEnergyDistribution{obj.N_Layer}.L(:,:,i,m+1);
                            T1(:,:,i,m+1) = obj.FullEnergyDistribution{i_layer}.T(:,:,i,m+1);
                            L1(:,:,i,m+1) = obj.FullEnergyDistribution{i_layer}.L(:,:,i,m+1);
                            Tw_cur(:,:,i,m+1) = obj.FullEnergyDistribution{i_layer}.T(:,:,i,m+1)*w + obj.FullEnergyDistribution{i_layer}.L(:,:,i,m+1);
                        end
                    end
                    A = conv3d(T1,wT2L, 'right',K_E)*obj.dE +conv3d(L1,T2L, 'right',K_E)*obj.dE;
                else
                    A = conv3d(Tw_cur,T,'right',K_E)*obj.dE;
                end
                C = conv3d(R,Rw_cur,'right',K_E)*obj.dE;
                T = obj.DoubleLayer(A,Tw_cur,C,T,K_E);
                %}
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
        
%         function CalculateAngularDistribution(obj,phi,theta0)
%             if nargin < 2; phi = []; end; %??????
%             if nargin < 3; theta0 = []; end;
%             obj.CalculateAngularDistribution_E0(phi,theta0);
%         end
        
    end % methods
end % classdef

