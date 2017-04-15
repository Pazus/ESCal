classdef NSAugerXPS < NSTransmition
    %NSXPS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Mat@XPSMaterial;
        Am
    end
    
    properties (SetAccess = protected, Transient)
%         F_m_m
%         F_p_m
        F_m
        F_p
        A
    end
    
    properties (Access = private, Transient)
        A_temp
        Ca_prev
        FAD
    end
    
    methods
        function obj = NSAugerXPS(Layer)
            obj = obj@NSTransmition(Layer);
            
            obj.CalculationResultPropertyName = 'Am';
        end
    end
    
    methods (Access = protected)
        
        function SolveForElasticScattering(obj)
            
            obj.SolveForElasticScattering@NSTransmition;
            
            obj.Ba = obj.AR0';
            C = obj.F_m + obj.F_p * obj.w * obj.R(:,:,1);
            
            if ~obj.isInfLayer
                obj.Ba = obj.Ba + obj.I;
                C = obj.h * C + obj.Ca_prev(:,:,1);
            end                    
                    
            obj.A(:,:,1) = C/obj.Ba;
        end
        
        function PrepareForInelasticScattegin(obj)
            obj.PrepareForInelasticScattegin@NSTransmition;
            
            obj.FAD = obj.F_p * obj.w + obj.A(:,:,1)*obj.D0;
        end
        
        function SolveForNInelasticScatterings(obj, i)
            obj.SolveForNInelasticScatterings@NSTransmition(i);

            Ca = obj.CalculateCa(i);
            obj.A(:,:,i) = Ca/obj.Ba;
        end
        
        function PreallocateResults(obj)
            obj.PreallocateResults@NSTransmition;
            obj.Am = zeros(obj.N, obj.N, obj.N_in+1, obj.M+1);
        end
        
        function PreallocateBuffers(obj)
            obj.PreallocateBuffers@NSTransmition;
            obj.A_temp = zeros(obj.N, obj.N, obj.N_in+1, obj.BDF_coefs.order);
            obj.A = zeros(obj.N, obj.N, obj.N_in+1);
            obj.FAD = zeros(obj.N);
        end
        
        function CalculateCrossSection(obj, varargin)
            obj.CalculateCrossSection@NSTransmition;

            norm = obj.Layer.Material.lambda * 1;
            obj.F_p = ones(obj.N)*norm/2/pi;
            obj.F_m = ones(obj.N)*norm/2/pi;
        end
        
        function CalculateOptimalM(obj, varargin)
            M = 1;% min(optimal_M(obj.theta0),2);
            if M ~= obj.M
                obj.M = M;
                obj.x_m_m = [];
                obj.x_p_m = [];
                obj.ClearResults;
            end;
        end
        
        function SaveResults(obj,m)
            obj.SaveResults@NSTransmition(m);
            obj.Am(:,:,:,m+1) = obj.A;
        end
        
        function Buffer(obj)
            obj.Buffer@NSTransmition;
            obj.A_temp = cat(4, obj.A, obj.A_temp(:,:,:,1:end-1));
        end
        
        function ClearPrevBuffer(obj)
            obj.Ca_prev = zeros(obj.N, obj.N, obj.N_in+1);
            obj.ClearPrevBuffer@NSTransmition;
        end
        
        function AppendToPrevBuffer(obj,BDF_coef,i)
            obj.Ca_prev = obj.Ca_prev + BDF_coef * obj.A_temp(:,:,:,i);
            obj.AppendToPrevBuffer@NSTransmition(BDF_coef,i);
        end
        
%         function BaseCalculation(obj)
%             lambda = obj.Layer.Material.lambda;
% 
%             for m = 0:0
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
%                     Ba = Ar';
%                     Ca = obj.F_m + obj.F_p * obj.w * obj.R(:,:,1);
%                     
%                     if ~obj.isInfLayer
%                         Ba = Ba + speye(obj.N)/2;
%                         Ca = obj.h * Ca + obj.Ca_prev(:,:,1);
%                     end
%                     
%                     obj.A(:,:,1) = Ca/Ba;
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
%                         obj.FAD = obj.F_p * obj.w + obj.A(:,:,1)*obj.D0;
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
%                             Ca = obj.CalculateCa;
%                             
%                             obj.A(:,:,i) = Ca/Ba;
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
%                 obj.Ca_prev = zeros(obj.N,obj.N,obj.N_in+1);
%                 % Сделал -1 для проверки! нужно обратить на это внимание!
%                 % В файлах Reflection2, Transmition2 для определения верхнего
%                 % предела цикла использовалось значение с -1, а в цикле бралось
%                 % без этой поправки. В данной реализации изменил, нужно
%                 % тестить.
%                 IterIntex = obj.TauIterationIndex;
%                 BDF_order = obj.BDF_coefs.order;
% 
%                 for i=1:min(BDF_order,IterIntex-1)
%                     obj.Ca_prev = obj.Ca_prev + obj.BDF_coefs.alpha(min(BDF_order,IterIntex),i) * obj.A_temp(:,:,:,i);
%                 end
%             end
%         end

    end
    
    methods (Access = private)
        function Ca = CalculateCa(obj, i)
            Ca = zeros(obj.N);
%             i = obj.InIterationIndex;
            x_in = obj.x_in;
            
            for j=2:floor(i/2)
                Ca = Ca + obj.A(:, :, j)*obj.DR(:, :, i-j);
            end
            Ca = Ca + Ca';
            
            if isodd(i)
                Ca = Ca + obj.A(:, :, ceil(i/2))*obj.DR(:, :, floor(i/2));
            end
            
            Ca = Ca + obj.A(:, :, i-1)*x_in;
            Ca = Ca + obj.FAD*obj.R(:, :, i);
            
            if ~obj.isInfLayer
                Ca = obj.h*Ca + obj.Ca_prev(:, :, i);
            end
        end
    end
    
end

