classdef PESMultiLayer<BaseMultiLayer
    properties
        %          z; %заменить все tau на z, вектор длиной в количество слоев
        theta                  = 0;
        phi                    = [];
        %         delta_max              = [];
        Fm;
        MatExist;
        ShellExist;
        Layer_flag;
        DIIMFP_Materials;
        calc_withoutR          = false;
    end
    
    methods
        function obj = PESMultiLayer (Layers, MatShells, CalculationMethod)
            
            
            % MatShells.MatName MatShells.ShellN0ame
            if nargin<3 || isempty(CalculationMethod);
                CalculationMethod = {'NS'};
            elseif nargin<2 %|| isempty(MatShells)
                %                 MatShells = % список элементов из Layers
                error('Set material shells')
            end
            
            
            % Проверка на одинаковость анодов
            % Проверка в слоях на наличие анализируемой оболочки и материала
            
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
            
            %{
            E_anode = Layers(1).Material.Anode.PhotonEnergy;
            MatData=fieldnames(MaterialData);
            if isempty(MatShells)
                MatShells = [];
                
            else
            %}
            MatExist=fieldnames(MatShells);
            obj.MatExist=MatExist;
            N_Mat = length(MatExist);
            
            for i_layer = 1:obj.N_Layer
                obj.DIIMFP_Materials(i_layer).DIIMFP_E = obj.Layers(i_layer).Material.DIIMFP_E;
                obj.DIIMFP_Materials(i_layer).DIIMFP = obj.Layers(i_layer).Material.DIIMFP;
                obj.DIIMFP_Materials(i_layer).E0 = obj.Layers(i_layer).Material.E0;
            end
            %{
                for i_mat = 1:length(MatExist)
                     if ~any(strcmp(MatExist,MatData));
                         error('Such material name is absent');
                     end
                end
            end
                
            
            
            if ~isfield(MatShells, 'ShellName')
                for i_mat = 1:length(MatShells.MatName)
                    FindMatProp=getfield(MaterialData, MatShells.MatName{i_mat});
                    obj.ShellsAnalyzed.MatName{i_mat} = MatShells.MatName{i_mat};
                    N_shells = length(FindMatProp.XPS.Shells);
                    k=1;
                    for i_sh = 1:N_shells
                        if  FindMatProp.XPS.EB(i_sh)<E_anode
                            obj.ShellsAnalyzed.ShellName{i_mat}(k)=FindMatProp.XPS.Shells(i_sh);
                            k=k+1;
                        end
                    end
                end
            else
                for i_mat = 1:length(MatShells.MatName)
                    FindMatProp=getfield(MaterialData, MatShells.MatName{i_mat});
                    if ~isempty(MatShells.ShellName{i_mat})
                        for i_shell=1:length(MatShells.ShellName{i_mat})
                            if any(strcmp(FindMatProp.XPS.Shells,MatShells.ShellName{i_mat}(i_shell)))
                                obj.ShellsAnalyzed.ShellName{i_mat}(i_shell)=MatShells.ShellName{i_mat}(i_shell);
                            else
                                error(['Неправильно заданы оболочки для ', MatShells.ShellName{i_mat} ])
                            end
                        end
                    else
                        obj.ShellsAnalyzed.MatName{i_mat} = MatShells.MatName{i_mat};
                        N_shells = length(FindMatProp.XPS.Shells);
                        k=1;
                        for i_sh = 1:N_shells
                            if  FindMatProp.XPS.EB<E_anode
                                obj.ShellsAnalyzed.ShellName{i_mat}(k)=FindMatProp.XPS.Shells(i_sh);
                                k=k+1;
                            end
                        end
                    end
                end
               
            end
            %}
            for i_N_Mat = 1:N_Mat
                obj.ShellExist.(MatExist{i_N_Mat}).Shell = MatShells.(MatExist{i_N_Mat});
                for j_Shell = 1:numel(obj.ShellExist.(MatExist{i_N_Mat}).Shell)
                    XPSMatLocal = XPSMaterial(MatExist{i_N_Mat},obj.ShellExist.(MatExist{i_N_Mat}).Shell{j_Shell});
                    obj.ShellExist.(MatExist{i_N_Mat}).E0(j_Shell) = XPSMatLocal.E0;
                end
                                
                obj.Layer_flag.(MatExist{i_N_Mat}) = false(1,obj.N_Layer);
                for i_layer = 1:obj.N_Layer
                    if min(Layers(i_layer).Material.DIIMFP_E)>5 || max(Layers(i_layer).Material.DIIMFP_E)<Layers(i_layer).Material.E0-5
                        error('DIIMFP should be set on interval 0:E0')
                    end
                    if strcmp(MatExist{i_N_Mat},Layers(i_layer).MaterialName)
                        obj.Layer_flag.(MatExist{i_N_Mat})(i_layer) = true;
                        
                    end
                end
            end
            
        end %PESMultiLayer
        
        function Calculate(obj)
            
            %             load MaterialData;
            for j_Mat=1:length(obj.MatExist)
                for j_Shell=1:length(obj.ShellExist.(obj.MatExist{j_Mat}).Shell)
                    %                     FindMatProp = getfield(MaterialData, obj.MatExist{j_Mat});
                    %                     Shell_number = find(strcmp(FindMatProp.XPS.Shells,obj.ShellExist.(obj.MatExist{j_Mat})(j_Shell)));
                    disp(obj.ShellExist.(obj.MatExist{j_Mat}).Shell{j_Shell})
                    for i_layer = 1:obj.N_Layer

                        if obj.Layer_flag.(obj.MatExist{j_Mat})(i_layer)

                            obj.Layers(i_layer).Material.SetShell(obj.ShellExist.(obj.MatExist{j_Mat}).Shell{j_Shell});
                            DIIMFP = obj.calcDIIMFP(i_layer,obj.Layers(i_layer).Material.E0);
                            obj.Layers(i_layer).Material.SetManualDIIMFP(obj.energy_mesh_full', DIIMFP');

                            switch obj.CalculationMethod{i_layer}
                                case 'NS'
                                    obj.ObjectsOfLayers{i_layer} = NSXPS(obj.Layers(i_layer));
                                case 'SLA'
                                    obj.ObjectsOfLayers{i_layer} = SLAXPS(obj.Layers(i_layer));
                                case 'SA'
                                    obj.ObjectsOfLayers{i_layer} = SAXPS(obj.Layers(i_layer));
                                otherwise
                                    error('Введите корректное название метода: NS, SLA или SA.');
                            end
                        else
%                             obj.Layers(find(obj.Layer_flag.(obj.MatExist{j_Mat}),true,'first')).Material.SetShell(obj.ShellExist.(obj.MatExist{j_Mat}){j_Shell});

                            Mat_local = Material(obj.Layers(i_layer).MaterialName,obj.ShellExist.(obj.MatExist{j_Mat}).E0(j_Shell));
                            DIIMFP = obj.calcDIIMFP(i_layer,obj.ShellExist.(obj.MatExist{j_Mat}).E0(j_Shell));
                            Mat_local.SetManualDIIMFP(obj.energy_mesh_full', DIIMFP');

                            Layer_local = Layer(Mat_local,obj.Layers(i_layer).thickness);
                            switch obj.CalculationMethod{i_layer}
                                case 'NS'
                                    obj.ObjectsOfLayers{i_layer} = NSTransmition(Layer_local);
                                case 'SLA'
                                    obj.ObjectsOfLayers{i_layer} = SLATransmition(Layer_local);
                                case 'SA'
                                    obj.ObjectsOfLayers{i_layer} = SATransmition(Layer_local);
                                otherwise
                                    error('Введите корректное название метода: NS, SLA или SA.');
                            end
                        end
                        
                    end
                    
                    obj.CalculateEachLayer;
                    obj.mu_mesh = obj.ObjectsOfLayers{1}.mu_mesh; % !!!!!
                    
                    obj.Setup_energy_mesh_full;
                    obj.RecalculateToEnergyDistribution;
%                     disp(['j_Mat=', num2str(j_Mat),', j_Shell])
                    if j_Mat == 1 && j_Shell == 1
                        obj.Fm = zeros(size(obj.FullEnergyDistribution{1}.T));
                    end
                    
                    if ~obj.Layer_flag.(obj.MatExist{j_Mat})(obj.N_Layer)
                        Q = [];
                    else
                        Q = obj.FullEnergyDistribution{obj.N_Layer}.Q;
                    end
                    
                    %             w = sparse(1:obj.ObjectsOfLayers{1}.N, 1:obj.ObjectsOfLayers{1}.N, obj.ObjectsOfLayers{1}.mu_mesh_weights./obj.mu_mesh);
                    if obj.N_Layer~=1
                        wT = zeros(size(obj.FullEnergyDistribution{1}.T));
                        Tw = zeros(size(obj.FullEnergyDistribution{1}.T));
                        R  = obj.FullEnergyDistribution{obj.N_Layer}.R;
                    end
                    M = size(obj.Fm,4)-1;
                    N_E =size(obj.Fm,3);
                    for i_layer = (obj.N_Layer-1):-1:1
                        if ~obj.Layer_flag.(obj.MatExist{j_Mat})(i_layer)
                            Qlocal = [];
                        else
                            Qlocal = obj.FullEnergyDistribution{i_layer}.Q(:,:,:,:);
                        end
                        %                 Rlocal = obj.FullEnergyDistribution{i_layer}.R(:,:,:,:);
                        if strcmp(obj.CalculationMethod{i_layer},'SLA')
                            % В SLA-решении уже содержатся частицы, прошедшие слой,
                            % не испытав упругих рассеяний, а сами решения
                            % представляют из себя дельта-функции, поэтому не нужно
                            % интегрировать и не нужно добавлять Решение Ландау,
                            % которое учитывает частицы, прошедшие слой, не испытав
                            % упругих рассеяний
                            wT = obj.FullEnergyDistribution{i_layer}.T;
                            Tw = obj.FullEnergyDistribution{i_layer}.T;
                            wR(:,:,:,:) = 1*R(:,:,:,:);
                        else
                            w = sparse(1:obj.ObjectsOfLayers{i_layer}.N, 1:obj.ObjectsOfLayers{i_layer}.N, obj.ObjectsOfLayers{i_layer}.mu_mesh_weights./obj.mu_mesh);
                            for m=0:M
                                for i=1:N_E
                                    wT(:,:,i,m+1) = w*obj.FullEnergyDistribution{i_layer}.T(:,:,i,m+1) + obj.FullEnergyDistribution{i_layer}.L(:,:,i,m+1);
                                    Tw(:,:,i,m+1) = obj.FullEnergyDistribution{i_layer}.T(:,:,i,m+1)*w + obj.FullEnergyDistribution{i_layer}.L(:,:,i,m+1);
                                    wR(:,:,i,m+1) = w*R(:,:,i,m+1);
                                end
                            end
                        end
                        E_temp=abs(obj.energy_mesh_full/obj.ObjectsOfLayers{i_layer}.Material.E0-1);
                        K_E = length(obj.energy_mesh_full)-find(min(E_temp)==E_temp,1); % Число элементов после E0
                        
                        if ~obj.calc_withoutR
                             Q_QR= obj.DoubleLayer(Q,     Qlocal,wR,[],K_E);
                             Q   = obj.DoubleLayer(Qlocal,Q_QR,  wT,[],K_E);
                        else Q   = obj.DoubleLayer(Qlocal,Q,     wT,[],K_E);
                        end
                        if i_layer<obj.N_Layer-1
                            Rlocal = obj.FullEnergyDistribution{i_layer}.R(:,:,:,:);
                            R = obj.DoubleLayer(Rlocal,Tw,R,wT,K_E);
                        end
                        
                    end
                    obj.Fm = obj.Fm + Q;
                end
            end
            obj.IsCalculated = true;
        end

        
        function y = calcDIIMFP(obj, i_layer, E0)
            dE=abs(obj.DIIMFP_Materials(i_layer).DIIMFP_E(2)-obj.DIIMFP_Materials(i_layer).DIIMFP_E(1));
            mesh_E = fliplr(E0:-dE:0);
            if E0>obj.DIIMFP_Materials(i_layer).E0
                delta_E0 = E0 - obj.DIIMFP_Materials(i_layer).E0;
                I_in_temp = interp1(obj.DIIMFP_Materials(i_layer).DIIMFP_E+delta_E0,obj.DIIMFP_Materials(i_layer).DIIMFP,mesh_E);
            elseif E0<obj.DIIMFP_Materials(i_layer).E0
                delta_E0 = obj.DIIMFP_Materials(i_layer).E0 - E0;
                I_in_temp = interp1(obj.DIIMFP_Materials(i_layer).DIIMFP_E-delta_E0,obj.DIIMFP_Materials(i_layer).DIIMFP,mesh_E);
            elseif E0==obj.DIIMFP_Materials(i_layer).E0
                delta_E0 = 0;
                I_in_temp = interp1(obj.DIIMFP_Materials(i_layer).DIIMFP_E-delta_E0,obj.DIIMFP_Materials(i_layer).DIIMFP,mesh_E);
            end
            I_in_temp(isnan(I_in_temp))=0;
            y = I_in_temp/trapz(mesh_E,I_in_temp);
            y = interp1(mesh_E,y,obj.energy_mesh_full);
            y(isnan(y))=0;
        end

        function CalculateEnergyDistribution(obj,theta,phi)
            if nargin < 3; phi = 0; end; %??????
            obj.CalculateEnergyDistribution_E0(theta,phi);
            convGauss (obj,obj.sigma_gauss);
            convLDS (obj,obj.sigma_LDS, obj.alpha_LDS);
        end
        
    end % methods
end % classdef