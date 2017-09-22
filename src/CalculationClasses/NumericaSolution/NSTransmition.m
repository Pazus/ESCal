classdef NSTransmition < NSReflection
    %NSTransmition Summary of this class goes here
    %   Detailed explanation goes here
    properties (SetAccess = protected)
        Tm
    end
    
    properties (Access = protected)
%         U
    end
    
    properties (Access = protected, Transient)
        current_tau
        exp_t0
        ETX
        CwR
        T
    end
    
    properties (Access = private, Transient)
        At
        Bt
        T_temp
        Ct_prev
    end
    
    methods
        function obj = NSTransmition(Layer)
            obj = obj@NSReflection(Layer);
            obj.CalculationResultPropertyName = 'Tm';
        end
    end
    
    methods (Access = protected)
        
        function SolveForElasticScattering(obj)
            lambda = obj.Layer.Material.lambda;
            
            obj.SolveForElasticScattering@NSReflection;
            
            if ~obj.isInfLayer
                obj.Bt = obj.AR0' + obj.I;
                obj.current_tau = obj.current_tau + obj.h0;

                % Косинуса рядом с экспонентой нет т.к. он сокращяется
                % Вероятность рассеятся в слое ~1/mu
                % Экспанента описывает поток, а плотность потока ~ mu
                % Косинусы сокращаются
                obj.exp_t0 = sparse(1:obj.N, 1:obj.N, exp(-obj.current_tau ./ obj.mu_mesh)); %./obj.x
                if obj.withoutNonliniearPart
                    C = obj.h * obj.exp_t0 * lambda * (obj.x_p) + obj.Ct_prev(:,:,1);
                else
                    C = obj.h * obj.exp_t0 * lambda * (obj.x_p + obj.x_m * obj.w * obj.R(:,:,1)) + obj.Ct_prev(:,:,1);
                end
                obj.T(:,:,1) = C/obj.Bt;
            end
        end
        
        function PrepareForInelasticScattegin(obj)
            obj.PrepareForInelasticScattegin@NSReflection;
            
            if ~obj.isInfLayer
                obj.ETX = (obj.exp_t0+obj.T(:,:,1)*obj.w)*obj.C0*obj.w;
            else
                obj.ETX = 0;
            end
        end
        function SolveForNInelasticScatterings(obj, i)
            obj.SolveForNInelasticScatterings@NSReflection(i);
            
            if obj.withoutNonliniearPart
                obj.CwR(:,:,i-1) = 0;
            else
                obj.CwR(:,:,i-1) = obj.C0*obj.w*obj.R(:,:,i);
            end
            
            if ~obj.isInfLayer
                Ct = obj.CalculateCt(i);

                obj.T(:,:,i) = Ct/obj.Bt;
            end;
        end
       
        function PreallocateResults(obj)
            obj.PreallocateResults@NSReflection;
            obj.Tm = zeros(obj.N,obj.N,obj.N_in+1,obj.M+1);
        end
        
        function PreallocateBuffers(obj)
            obj.PreallocateBuffers@NSReflection;
            obj.T_temp = zeros(obj.N,obj.N,obj.N_in+1,obj.BDF_coefs.order);
            obj.T = zeros(obj.N,obj.N,obj.N_in+1);
            obj.CwR = zeros(obj.N,obj.N,obj.N_in);
        end
        
        function SaveResults(obj,m)
            obj.SaveResults@NSReflection(m);
            obj.Tm(:,:,:,m+1) = obj.T;
        end
        
        function Buffer(obj)
            obj.Buffer@NSReflection;
            obj.T_temp = cat(4,obj.T,obj.T_temp(:,:,:,1:end-1));
        end
        
        function ClearPrevBuffer(obj)
            obj.Ct_prev = zeros(obj.N,obj.N,obj.N_in+1);
            obj.ClearPrevBuffer@NSReflection;
        end
        function AppendToPrevBuffer(obj,BDF_coef,i)
            obj.Ct_prev = obj.Ct_prev + BDF_coef * obj.T_temp(:,:,:,i);
            obj.AppendToPrevBuffer@NSReflection(BDF_coef,i);
        end

        function BasePrecalculation(obj)
            obj.BasePrecalculation@NSReflection;
            obj.current_tau = 0;
        end
        
%         function BaseCalculation(obj)
%             
%             if obj.isInfLayer
%                 obj.BaseCalculation@NSReflection;
%                 return;
%             end
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
%                     Cr = obj.h * obj.C0 + obj.Cr_prev(:,:,1);
%                     Ar = obj.h * (diag(1./obj.mu_mesh) - obj.G) + obj.I;
%                     Dr = obj.h * obj.D0;
%                     
%                     obj.R(:,:,1) = symmetrize(care_Newton(Ar, Cr, Dr, obj.eps, obj.max_newton_iter_count));
%                     
%                     Ar = Ar - obj.R(:,:,1)*Dr;
%                     
%                     obj.current_tau = obj.h0 * (k-1);
%                     % Косинуса рядом с экспонентой нет т.к. он сокращяется
%                     % Вероятность рассеятся в слое ~1/mu
%                     % Экспанента описывает поток, а плотность потока ~ mu
%                     % Косинусы сокращаются
%                     obj.exp_t0 = sparse(1:obj.N, 1:obj.N, exp(-obj.current_tau ./ obj.mu_mesh)); %./obj.x
%                     
% %                     Bt = obj.h*(diag(1 ./ obj.mu_mesh) - obj.G' - obj.D0*obj.R(:,:,1)) + speye(obj.N); % diag(1./x)
%                     Bt = Ar' + obj.I; % diag(1./x)
%                     Ct = obj.h * obj.exp_t0 * lambda * (obj.x_p + obj.x_m * obj.w * obj.R(:,:,1)) + obj.Ct_prev(:,:,1);
%                     
%                     obj.T(:,:,1) = Ct/Bt;
% 
%                     if obj.N_in > 0
%                         % Пересчитываем для кратности неупругого рассеяния большей нуля
%                         
%                         % Находим разложение по собственным числам и собственным
%                         % векторам постоянных коэффициентов уравнения
%                         [Ur,Ar_eig] = eig_sorted(Ar);
%                         obj.ETX = (obj.exp_t0+obj.T(:,:,1)*obj.w)*obj.C0*obj.w;
%                         
%                         
%                         for i=2:obj.N_in+1
%                             obj.InIterationIndex = i;
%                             
%                             Cr = obj.CalculateCr;
%                             
%                             obj.R(:,:,i) = sylvsolve_eig(Ar_eig, Cr, Ur);
%                             obj.DR(:,:,i-1) = obj.D0 * obj.R(:,:,i);
%                             obj.CwR(:,:,i-1) = obj.C0*obj.w*obj.R(:,:,i);
%                             
%                             Ct = obj.CalculateCt;
%                             
%                             obj.T(:,:,i) = Ct/Bt;
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
%                 obj.CalculatePrevBuffers@NSReflection;
%                 obj.Ct_prev = zeros(obj.N,obj.N,obj.N_in+1);
%                 % Сделал -1 для проверки! нужно обратить на это внимание!
%                 % В файлах Reflection2, Transmition2 для определения верхнего
%                 % предела цикла использовалось значение с -1, а в цикле бралось
%                 % без этой поправки. В данной реализации изменил, нужно
%                 % тестить.
%                 IterIntex = obj.TauIterationIndex;
%                 BDF_order = obj.BDF_coefs.order;
% 
%                 for i=1:min(BDF_order,IterIntex-1)
%                     obj.Ct_prev = obj.Ct_prev + obj.BDF_coefs.alpha(min(BDF_order,IterIntex),i) * obj.T_temp(:,:,:,i);
%                 end
%             end
%         end

        
    end
    
    methods (Access = private) 
        function Ct = CalculateCt(obj, i)

            Ct = zeros(obj.N);
%             i = obj.InIterationIndex;
            w = obj.w;
            x_in = obj.x_in;

            CwR_inv = flip(obj.CwR(:,:,1:i-2),3);

            exp_t = obj.exp_t0;
            for j=2:i-1
                exp_t = exp_t.*(obj.current_tau*x_in/(j-1));
                Ct = Ct + (exp_t+obj.T(:,:,j)*w)*CwR_inv(:,:,j-1);
            end
            
            
            if obj.withoutNonliniearPart
                Ct = Ct + (exp_t.*diag(obj.mu_mesh) + obj.T(:,:,i-1))*x_in  ... 
                    + exp_t.*(obj.current_tau*x_in/(i-1))*obj.Layer.Material.lambda*(obj.x_p);
            else
                Ct = Ct + (exp_t.*diag(obj.mu_mesh) + obj.T(:,:,i-1))*x_in  ... 
                    + exp_t.*(obj.current_tau*x_in/(i-1))*obj.Layer.Material.lambda*(obj.x_p+obj.x_m*w*obj.R(:,:,1)) ...
                    + obj.ETX * obj.R(:,:,i);
            end

            Ct = obj.h*Ct + obj.Ct_prev(:,:,i);

        end
    end
    
end

