classdef NSBase < BaseLayerCalculation
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = protected, SetObservable = true)
%         TauIterationIndex;                      % Current iteration for tau
        Nt                                      % numel tau mesh
    end
    
    properties (SetAccess = protected)
        w                                       % diag(s./x)

%         InIterationIndex;                       % Current iteration for N_ins
%         MIterationIndex;                        % Current iteration for M
        BDF_coefs@BDF                           % Coefficients of BDF method
        
        x_m_m                                   % Recalculated transmition DECS for all m values
        x_p_m                                   % Recalculated reflection DECS for all m values
        x_m                                     % Recalculated transmition DECS for current m value
        x_p                                     % Recalculated reflection DECS for current m value
        
        IsManualNt              = false;        % Manual Nt value assigned flag
    end
    
    properties (Access = protected)
        h0                                      % tau mesh step
        h                                       % BDF method step
        isInfLayer = true;
    end
    
    properties (Constant = true)
        r                       = 6;            % order of BDF method
        eps                     = 1e-10;        % Riccati equation solution error
        max_newton_iter_count   = 100;          % number of Newton iterations before failure
    end
    
    properties (Dependent = true, SetAccess = protected)
        I                                       % Precalculation coef I = eye(N)/2
    end
    
    methods
        function obj = NSBase(Layer)
            
            obj = obj@BaseLayerCalculation(Layer);
            
%             addlistener(obj,'TauIterationIndex','PostSet',@obj.CalculateH);
            addlistener(obj,'Nt','PostSet',@obj.ClearResults);
            addlistener(obj.Layer,'LayerRecalculated',@obj.CalculateOptimalNt);
            
            obj.BDF_coefs = BDF(obj.r);
            
            obj.CalculateOptimalNt;
%             obj.CalculateCrossSection;
        end
        
        function CalculateOptimalNt(obj, varargin)
            if ~obj.IsManualNt || isempty(varargin)
                if isinf(obj.Layer)
                    obj.Nt = 1;
                else
                    obj.Nt = ceil(min(max(obj.Layer.tau_tot*12,15),100));
                end
                obj.IsManualNt = false;
            end
        end
        
        function SetManualNt(obj,Nt)
            ps = inputParser;
            ps.FunctionName = 'SetManualNt';
            ps.addRequired('Nt',@(x) validateattributes(x,{'numeric'},{'integer','nonnegative','scalar','nonempty'}));
            ps.parse(Nt);
            
            obj.Nt = ps.Results.Nt;
            obj.IsManualNt = true;
        end
        
        function Calculate(obj)
            obj.ClearResults;
            
            obj.CalculateAngleMesh;
            obj.CalculateCrossSection;
            
            obj.PreallocateResults;
            
            obj.isInfLayer = isinf(obj.Layer);
            
            obj.BaseCalculation2;
            
            obj.IsCalculated = true;
        end
        
        function BaseCalculation2(obj)
            obj.CalculateH0;
            for m=0:obj.M
                
                obj.SetCrossection(m);
                obj.PreallocateBuffers;
                obj.BasePrecalculation;
                
                for k = 1:obj.Nt
                    % Solution for N0
                    
%                     obj.TauIterationIndex = k;
%                     obj.InIterationIndex = 1;
                    obj.CalculateH(k);
                    
                    obj.CalculatePrevBuffers2(k);
                    
                    obj.SolveForElasticScattering;
                    
                    if obj.N_in > 0
                        % Fill buffering properties
                        obj.PrepareForInelasticScattegin;
                        
                        for i=2:obj.N_in+1
%                             obj.InIterationIndex = i;
                            
                            obj.SolveForNInelasticScatterings(i);
                        end
                    end
                    
                    obj.Buffer;
                end
                
                obj.SaveResults(m);
                
            end
        end
        
        function val = get.I(obj)
            if ~obj.isInfLayer
                val = speye(obj.N)/2;
            else
                val = 0;
            end
        end

    end
    methods (Access = protected)
        
        function ClearResults(obj,varargin)
            obj.x_m        = [];
            obj.x_p        = [];
            
            obj.ClearResults@BaseCalculation;
        end
        
        function CalculateCrossSection(obj,varargin)
            
%             [obj.x_p_m, obj.x_m_m] = X_mp(obj.Layer.Material.DECS,obj.Layer.Material.DECS_mu,obj.mu_mesh,obj.mu_mesh_weights,obj.M,1);
            x_l = obj.Layer.Material.Calculate_Leg_coefs(2000)*obj.norm;
            [obj.x_p_m, obj.x_m_m] = expandLegCoefsTo2ParamCrossSection(x_l, obj.mu_mesh,obj.M);
        end
        
        function CalculateH0(obj)
            if obj.isInfLayer
                obj.h0 = 1;
            else
                obj.h0 = obj.Layer.tau_tot / (obj.Nt);
            end
        end
        
        function CalculateH(obj, k)
            if obj.isInfLayer
                obj.h = 1;
            else
                obj.h = obj.h0 * obj.BDF_coefs.betta(min(k,obj.BDF_coefs.order));
            end
        end
        
        function CalculateAngleMesh(obj, varargin)
            obj.CalculateAngleMesh@BaseCalculation;
            obj.w = sparse(1:obj.N,1:obj.N,obj.mu_mesh_weights./obj.mu_mesh);
        end
        
        function CalculatePrevBuffers2(obj,k)
            if ~obj.isInfLayer
                obj.ClearPrevBuffer;
                
                % —делал -1 дл€ проверки! нужно обратить на это внимание!
                % ¬ файлах Reflection2, Transmition2 дл€ определени€ верхнего
                % предела цикла использовалось значение с -1, а в цикле бралось
                % без этой поправки. ¬ данной реализации изменил, нужно
                % тестить.
                BDF_order = obj.BDF_coefs.order;

                for i=1:min(BDF_order,k-1)
                    BDF_coef = obj.BDF_coefs.alpha(min(BDF_order,k),i);
                    obj.AppendToPrevBuffer(BDF_coef,i);
                end
            end
        end
        
    end
    
    methods (Access = protected, Abstract)
        PreallocateResults(obj);
%         BaseCalculation(obj);
        ClearPrevBuffer(obj);
        BasePrecalculation(obj);
        PreallocateBuffers(obj);
        SetCrossection(obj, m);
        AppendToPrevBuffer(obj, k, i);
        SolveForElasticScattering(obj);
        PrepareForInelasticScattegin(obj);
        Buffer(obj);
        SaveResults(obj, m);
        
        SolveForNInelasticScatterings(obj, i);
    end
end