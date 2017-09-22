classdef NSReflection < NSBase
    %NSReflection Solution of elastic electron reflection
    %   Equation is solved for R
    %   R_temp is a buffer for solution of d < Inf
    %
    %   Equation form for N_in = 0:
    %
    %       A*R+R*A' = C0 + R*D*R
    %
    %   Equation form for N_in = 0:
    %
    %       A*R+R*A' = Cn
    %
    properties (SetAccess = protected)
        Rm                                       % Result
    end
    
    
    properties (Access = protected, Transient)
        %% Precalculation coefs
        G                                       % Precalculation coef G
        C0                                      % Precalculation coef C0
        D0                                      % Precalculation coef D0
        x_in                                    % Precalculation coef x_in
        %% Buffer
        DR
        R
        AR0
        %% plotting config
    end
    
    properties (Access = private, Transient)
        U_eig
        A_eig
        Cr_prev
        R_temp                                  % Buffer of solutions. Used on d < Inf
    end
    
    properties 
        withoutNonliniearPart = false;
    end
    
    methods

        function obj = NSReflection(Layer)
            obj = obj@NSBase(Layer);
            obj.CalculationResultPropertyName   = 'Rm';
        end

    end
    methods (Access = protected)
        
        function SolveForElasticScattering(obj)
            if ~obj.isInfLayer
                Cr = obj.h * obj.C0 + obj.Cr_prev(:,:,1);
                Ar = obj.h * (diag(1./obj.mu_mesh) - obj.G) + obj.I;
                Dr = obj.h * obj.D0;
            else
                Cr = obj.C0;
                Ar = diag(1./obj.mu_mesh) - obj.G;
                Dr = obj.D0;
            end

            obj.R(:,:,1) = symmetrize(care_Newton(Ar,Cr,Dr,obj.eps,obj.max_newton_iter_count));
%             [Ud,D_eig] = eig_sorted(Dr);
%             obj.R(:,:,1) = symmetrize(care(Ar,Ud,Cr,diag(D_eig.^-1),zeros(obj.N),eye(obj.N)));
            
            obj.AR0 = Ar - obj.R(:,:,1)*Dr;
        end
        
        function PrepareForInelasticScattegin(obj)
            [obj.U_eig, obj.A_eig] = eig_sorted(obj.AR0);
        end
        
        function SolveForNInelasticScatterings(obj, i)
            Cr = obj.CalculateCr(i);

            obj.R(:,:,i) = sylvsolve_eig(obj.A_eig, Cr, obj.U_eig);
            obj.DR(:,:,i-1) = obj.D0 * obj.R(:,:,i);
        end
        
        function BasePrecalculation(obj)
            lambda = obj.Layer.Material.lambda;
            w = obj.w;
            x = obj.mu_mesh;
            N = obj.N;
            
            obj.G = lambda * obj.x_p * w;
            obj.x_in = sparse(1:N, 1:N, (1-lambda) ./ x);
            
            obj.C0 = lambda * obj.x_m;
            if obj.withoutNonliniearPart
                obj.D0 = 0;
            else
                obj.D0 = symmetrize(w * obj.C0 * w);
            end
%             obj.D0 = symmetrize(w * obj.C0 * w);
        end
        
        function SetCrossection(obj,m)
            obj.x_m = obj.x_m_m(:,:,m+1);
            obj.x_p = obj.x_p_m(:,:,m+1);
        end
        
        function PreallocateResults(obj)
            obj.Rm = zeros(obj.N,obj.N,obj.N_in+1,obj.M+1);
        end
        
        function PreallocateBuffers(obj)
            obj.R_temp = zeros(obj.N,obj.N,obj.N_in+1,obj.BDF_coefs.order);
            obj.DR = zeros(obj.N,obj.N,obj.N_in);
            obj.R = zeros(obj.N,obj.N,obj.N_in+1);
        end
        
        function Buffer(obj)
            if ~obj.isInfLayer
                obj.R_temp = cat(4,obj.R,obj.R_temp(:,:,:,1:end-1));
            end
        end
        
        function SaveResults(obj,m)
            obj.Rm(:,:,:,m+1) = obj.R;
        end
        
        function ClearResults(obj,varargin)
            obj.ClearResults@NSBase;
            obj.G          = [];
            obj.C0         = [];
            obj.D0         = [];
            obj.x_in       = [];
            obj.R_temp     = [];
            obj.DR         = [];
            obj.Cr_prev    = [];
            obj.R          = [];
            obj.U_eig      = [];
            obj.A_eig      = [];
            obj.Rm         = [];
        end
        
        function ClearPrevBuffer(obj)
            obj.Cr_prev = zeros(obj.N,obj.N,obj.N_in+1);
        end
        function AppendToPrevBuffer(obj,BDF_coef,i)
            obj.Cr_prev = obj.Cr_prev + BDF_coef * obj.R_temp(:,:,:,i);
        end
            
%         function BaseCalculation(obj)
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
%                         Cr = obj.h * obj.C0 + obj.Cr_prev(:,:,1);
%                         Ar = obj.h * (diag(1./obj.mu_mesh) - obj.G) + obj.I;
%                         Dr = obj.h * obj.D0;
%                     else
%                         Cr = obj.C0;
%                         Ar = diag(1./obj.mu_mesh) - obj.G;
%                         Dr = obj.D0;
%                     end
%                     
%                     obj.R(:,:,1) = symmetrize(care_Newton(Ar,Cr,Dr,obj.eps,obj.max_newton_iter_count));
%                     
%                     if obj.N_in > 0
%                         % Пересчитываем для кратности неупругого рассеяния большей нуля
%                         Ar = Ar - obj.R(:,:,1)*Dr;
%                         % Находим разложение по собственным числам и собственным
%                         % векторам постоянных коэффициентов уравнения
%                         [Ur,Ar_eig] = eig_sorted(Ar);
%                         
%                         for i=2:obj.N_in+1
%                             obj.InIterationIndex = i;
%                             
%                             Cr = obj.CalculateCr;
%                             
%                             obj.R(:,:,i) = sylvsolve_eig(Ar_eig, Cr, Ur);
%                             obj.DR(:,:,i-1) = obj.D0 * obj.R(:,:,i);
%                         end
%                         
%                     end
%                     obj.Buffer;
%                 end
%                 obj.SaveResults(m);
%             end
%             
%         end
%         
%         function SolveEquation(obj)
%             n_in = obj.InIterationIndex;
%             if n_in == 1
%                 obj.R(:,:,n_in) = care_Newton(obj.Ar,obj.Cr,obj.Dr,obj.eps,obj.max_newton_iter_count);
%             else
%                 obj.R(:,:,n_in) = sylvsolve_eig(obj.Ar_eig,obj.Cr,obj.Ur);
%                 obj.DR(:,:,n_in-1) = obj.Dr * obj.R(:,:,n_in);
%             end
%         end
%
%         function Precalculation2(obj)
%             % Пересчитываем для кратности неупругого рассеяния большей нуля
%             obj.Ar = obj.Ar - obj.R(:,:,1)*obj.Dr;
%             % Находим разложение по собственным числам и собственным
%             % векторам постоянных коэффициентов уравнения
%             [obj.Ur,obj.Ar_eig] = eig_sorted(obj.Ar);
%         end
%
%         function CalculatePrevBuffers(obj)
%             if ~obj.isInfLayer
%                 obj.Cr_prev = zeros(obj.N,obj.N,obj.N_in+1);
%                 % Сделал -1 для проверки! нужно обратить на это внимание!
%                 % В файлах Reflection2, Transmition2 для определения верхнего
%                 % предела цикла использовалось значение с -1, а в цикле бралось
%                 % без этой поправки. В данной реализации изменил, нужно
%                 % тестить.
%                 IterIntex = obj.TauIterationIndex;
%                 BDF_order = obj.BDF_coefs.order;
% 
%                 for i=1:min(BDF_order,IterIntex-1)
%                     obj.Cr_prev = obj.Cr_prev + obj.BDF_coefs.alpha(min(BDF_order,IterIntex),i) * obj.R_temp(:,:,:,i);
%                 end
%             end
%         end
%
%         function CalculateEquationCoefs(obj)
%             if ~obj.isInfLayer
%                 h = obj.h;%_current;
%                 obj.Cr = h * obj.C0 + obj.Cr_prev(:,:,1);
%                 obj.Ar = h * (diag(1./obj.mu_mesh) - obj.G) + speye(obj.N)/2;
%                 obj.Dr = h * obj.D0;
%             else
%                 % Пользуемся преимущество Layzy copy.
%                 % Т.к. мы просто присваиваем значение другой переменной,
%                 % физически ничего в памяти не копируется, а просто
%                 % ставится поинтер.
%                 obj.Cr = obj.C0;
%                 obj.Ar = diag(1./obj.mu_mesh) - obj.G;
%                 obj.Dr = obj.D0;
%             end
%         end
        
    end
    
    methods (Access = private)
        function Cr = CalculateCr(obj, i)
%             i = obj.InIterationIndex;
            Cr = zeros(obj.N);

            for jj=2:floor(i/2)
                Cr = Cr + obj.R(:,:,jj) * obj.DR(:,:,i-jj);
            end
            Cr = Cr + Cr';

            if mod(i, 2) == 1 % isodd(i) perfirmance gain
                Cr = Cr + obj.R(:,:,ceil(i/2)) * obj.DR(:,:,floor(i/2));
            end

            x_in_R = obj.x_in * obj.R(:,:,i-1);
            Cr = Cr + x_in_R+x_in_R';

            if ~obj.isInfLayer
                Cr = obj.h * Cr + obj.Cr_prev(:,:,i);
            end
        end
    end
end

