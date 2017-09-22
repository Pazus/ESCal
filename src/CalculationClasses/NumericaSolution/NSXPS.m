classdef NSXPS < NSTransmition
    %NSXPS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Mat@XPSMaterial;
        Qm
        flipLayer = false;
    end
    
    properties (SetAccess = protected, Transient)
        F_m_m
        F_p_m
        F_m
        F_p
        Q
    end
    
    properties (Access = private, Transient)
        Cq_prev
        FQD
        Bq
        Q_temp
    end
    
    methods
        function obj = NSXPS(Layer)
            if ~isa(Layer.Material,'XPSMaterial'); error('Material of the Layer variable must be an "XPSMaterial" class object'); end
            obj = obj@NSTransmition(Layer);
            
            obj.CalculationResultPropertyName = 'Qm';
        end
    end
    
    methods (Access = protected)
        function SolveForElasticScattering(obj)
            
            obj.SolveForElasticScattering@NSTransmition;
            
            obj.Bq = obj.AR0';
            if obj.withoutNonliniearPart
                C = obj.F_m;
            else
                C = obj.F_m + obj.F_p * obj.w * obj.R(:,:,1);
            end
            
            if ~obj.isInfLayer
                obj.Bq = obj.Bq + obj.I;
                C = obj.h * C + obj.Cq_prev(:,:,1);
            end                    
                    
            obj.Q(:,:,1) = C/obj.Bq;
        end
        
        function PrepareForInelasticScattegin(obj)
            obj.PrepareForInelasticScattegin@NSTransmition;
            
            obj.FQD = obj.F_p * obj.w + obj.Q(:,:,1)*obj.D0;
        end
        
        function SolveForNInelasticScatterings(obj, i)
            obj.SolveForNInelasticScatterings@NSTransmition(i);

            Cq = obj.CalculateCq(i);
            obj.Q(:,:,i) = Cq/obj.Bq;
        end
        
        function PreallocateResults(obj)
            obj.PreallocateResults@NSTransmition;
            obj.Qm = zeros(obj.N, obj.N, obj.N_in+1, obj.M+1);
        end
        
        function PreallocateBuffers(obj)
            obj.PreallocateBuffers@NSTransmition;
            obj.Q_temp = zeros(obj.N, obj.N, obj.N_in+1, obj.BDF_coefs.order);
            obj.Q = zeros(obj.N, obj.N, obj.N_in+1);
            obj.FQD = zeros(obj.N);
        end
        
        function CalculateCrossSection(obj, varargin)
            obj.CalculateCrossSection@NSTransmition;
            x = -1:0.001:1;
            
            F = obj.Layer.Material.CalculateDPCS(x);
            
            if obj.flipLayer
                F = flip(F);
            end
%             P = P_Leg(x,3);
%             F = 1 - obj.Layer.Material.betta/2*P(:,3)';
            norm = obj.Layer.Material.lambda*obj.Layer.Material.lambda_photon;%obj.norm*
            [obj.F_p_m, obj.F_m_m] = X_mp(F, x, obj.mu_mesh, obj.mu_mesh_weights, obj.M,norm);
        end
        
%         function CalculateOptimalM(obj, varargin)
%             M = min(optimal_M(obj.theta0),2);
%             if M ~= obj.M
%                 obj.M = M;
%                 obj.x_m_m = [];
%                 obj.x_p_m = [];
%                 obj.ClearResults;
%             end;
%         end

        function SetCrossection(obj, m)
            obj.SetCrossection@NSTransmition(m);
            obj.F_m = obj.F_m_m(:,:,m+1);
            obj.F_p = obj.F_p_m(:,:,m+1);
        end
        
        function SaveResults(obj,m)
            obj.SaveResults@NSTransmition(m);
            obj.Qm(:,:,:,m+1) = obj.Q;
        end
        
        function Buffer(obj)
            obj.Buffer@NSTransmition;
            obj.Q_temp = cat(4, obj.Q, obj.Q_temp(:,:,:,1:end-1));
        end
        
        function ClearPrevBuffer(obj)
            obj.Cq_prev = zeros(obj.N, obj.N, obj.N_in+1);
            obj.ClearPrevBuffer@NSTransmition;
        end
        
        function AppendToPrevBuffer(obj,BDF_coef,i)
            obj.Cq_prev = obj.Cq_prev + BDF_coef * obj.Q_temp(:,:,:,i);
            obj.AppendToPrevBuffer@NSTransmition(BDF_coef,i);
        end
        
%         function BaseCalculation(obj)
%             
%             lambda = obj.Layer.Material.lambda;
% 
%             for m = 0:obj.M
% %                 obj.MIterationIndex = m;
%                 
%                 obj.SetCrossection(m);
%                 obj.PreallocateBuffers;
%                 obj.BasePrecalculation;
%                 
%                 for k = 1:obj.Nt
%                     obj.TauIterationIndex = k;
%                     obj.InIterationIndex = 1;
%                     
%                     obj.CalculatePrevBuffers;
%                     
%                     if ~obj.isInfLayer
%                         Ar = obj.h * (diag(1./obj.mu_mesh) - obj.G) + speye(obj.N)/2;
%                         Cr = obj.h * obj.C0 + obj.Cr_prev(:,:,1);
%                         Dr = obj.h * obj.D0;
%                     else
%                         Ar = diag(1./obj.mu_mesh) - obj.G;
%                         Cr = obj.C0;
%                         Dr = obj.D0;
%                     end
%                     
%                     obj.R(:,:,1) = symmetrize(care_Newton(Ar,Cr,Dr,obj.eps,obj.max_newton_iter_count));
%                     
%                     Ar = Ar - obj.R(:,:,1)*Dr;
%                     
%                     Bq = Ar';
%                     Cq = obj.F_m + obj.F_p * obj.w * obj.R(:,:,1);
%                     
%                     if ~obj.isInfLayer
%                         Bq = Bq + speye(obj.N)/2;
%                         Cq = obj.h * Cq + obj.Cq_prev(:,:,1);
%                     end
%                     
%                     obj.Q(:,:,1) = Cq/Bq;
%                     
%                     if ~obj.isInfLayer
%                         obj.current_tau = obj.h0 * (k);
%                         obj.exp_t0 = sparse(1:obj.N, 1:obj.N, exp(-obj.current_tau ./ obj.mu_mesh)); %./obj.x
% 
%     %                     Bt = obj.h*(diag(1 ./ obj.x) - obj.G' - obj.D0*obj.R(:,:,1)) + speye(obj.N); % diag(1./x)
%                         Bt = Ar' + speye(obj.N)/2; % diag(1./x)
%                         Ct = obj.h * obj.exp_t0 * lambda * (obj.x_p + obj.x_m * obj.w * obj.R(:,:,1)) + obj.Ct_prev(:,:,1);
% 
%                         obj.T(:,:,1) = Ct/Bt;
%                     end
%                     
%                     if obj.N_in > 0
%                         % Пересчитываем для кратности неупругого рассеяния большей нуля
%                         
%                         % Находим разложение по собственным числам и собственным
%                         % векторам постоянных коэффициентов уравнения
%                         [Ur,Ar_eig] = eig_sorted(Ar);
%                         obj.FQD = obj.F_p * obj.w + obj.Q(:,:,1)*obj.D0;
%                         
%                         if ~obj.isInfLayer
%                             obj.ETX = (obj.exp_t0+obj.T(:,:,1)*obj.w)*obj.C0*obj.w;
%                         end
%                         
%                         for i=2:obj.N_in+1
%                             obj.InIterationIndex = i;
%                             
%                             Cr = obj.CalculateCr;
%                             
%                             obj.R(:,:,i) = sylvsolve_eig(Ar_eig, Cr, Ur);
%                             obj.DR(:,:,i-1) = obj.D0 * obj.R(:,:,i);
%                             
%                             Cq = obj.CalculateCq;
%                             
%                             obj.Q(:,:,i) = Cq/Bq;
%                             
%                             if ~obj.isInfLayer
%                             
%                                 obj.CwR(:,:,i-1) = obj.C0*obj.w*obj.R(:,:,i);
%                                 Ct = obj.CalculateCt;
%                                 
%                                 obj.T(:,:,i) = Ct/Bt;
%                             end
%                         end
%                         
%                     end
%                     obj.Buffer;
%                 end
%                 obj.SaveResults(m);
%             end
%         end
% 
%         function CalculatePrevBuffers(obj)
%             if ~obj.isInfLayer
%                 obj.CalculatePrevBuffers@NSTransmition;
%                 obj.Cq_prev = zeros(obj.N,obj.N,obj.N_in+1);
%                 % Сделал -1 для проверки! нужно обратить на это внимание!
%                 % В файлах Reflection2, Transmition2 для определения верхнего
%                 % предела цикла использовалось значение с -1, а в цикле бралось
%                 % без этой поправки. В данной реализации изменил, нужно
%                 % тестить.
%                 IterIntex = obj.TauIterationIndex;
%                 BDF_order = obj.BDF_coefs.order;
% 
%                 for i=1:min(BDF_order,IterIntex-1)
%                     obj.Cq_prev = obj.Cq_prev + obj.BDF_coefs.alpha(min(BDF_order,IterIntex),i) * obj.Q_temp(:,:,:,i);
%                 end
%             end
%         end


    end
    
    methods (Access = private)
        function Cq = CalculateCq(obj, i)
            Cq = zeros(obj.N);
%             i = obj.InIterationIndex;
            x_in = obj.x_in;
            
            for j=2:floor(i/2)
                Cq = Cq + obj.Q(:,:,j)*obj.DR(:,:,i-j);
            end
            Cq = Cq + Cq';
            
            if isodd(i)
                Cq = Cq + obj.Q(:,:,ceil(i/2))*obj.DR(:,:,floor(i/2));
            end
            
            Cq = Cq + obj.Q(:,:,i-1)*x_in;
            Cq = Cq + obj.FQD*obj.R(:,:,i);
            
            if ~obj.isInfLayer
                Cq = obj.h*Cq + obj.Cq_prev(:,:,i);
            end
        end
    end
    
end

