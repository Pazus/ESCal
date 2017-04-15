classdef SABase < BaseLayerCalculation
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    
    properties (SetObservable = true, AbortSet = true)
        LegandePolinomialsCount   = 100;            % number of legandre polinomials for calculation
        CSCalculationMeshSizeCoef  = 5;
    end
    
    properties (SetAccess = protected)
        mu_full
        s_full
        x_l
    end
    
    properties (Dependent = true, SetAccess = private)
        N_full
    end
    
    properties (Access = protected, Transient)
        mu_mesh_temp
        mu_weight_temp
    end
    
    methods
        function obj = SABase(Layer)
            obj = obj@BaseLayerCalculation(Layer);
            
            addlistener(obj,'LegandePolinomialsCount','PostSet',@obj.ClearResults);
            addlistener(obj,'CSCalculationMeshSizeCoef','PostSet',@obj.ClearResults);
            
        end
        
        function Calculate(obj)
            obj.ClearResults;
            obj.CalculateCrossSection;
            
            obj.BaseCalculation;
            
            obj.AllocateResults;
            
            obj.SummarizeExpansionCoefficients;
            
            obj.IsCalculated = 1;
        end
        
        function val = get.N_full(obj)
            val = 2*obj.N-1;
        end

    end
    
    methods (Access = protected)
        function CalculateAngleMesh(obj,varargin)
            obj.CalculateAngleMesh@BaseCalculation;
            
            [x, s] = legzo_n1(obj.N_full-2);
            obj.mu_full = [1, x, -1]; obj.s_full = [0, s, 0];
        end;
        function CalculateCrossSection(obj, varargin)
            
            obj.CalculateTempMeshForExpansion;
            obj.ExpandCrossSections;
            
        end
        
        function SummarizeExpansionCoefficients(obj)
           norm_leg = ((0:obj.LegandePolinomialsCount)*2+1)/2;
           for m=0:obj.M
                P = Legendre_mu(obj.mu_mesh, m, obj.LegandePolinomialsCount)';
                for i=1:obj.LegandePolinomialsCount+1-m
                    PlP = P(:,i)*norm_leg(i+m)*P(:,i)'; %
                    for n=0:obj.N_in
                        obj.AddExpansionCoef(n,m,i,PlP);
                    end
                end
            end
        end
    end
    
    methods (Access = private)
        function CalculateTempMeshForExpansion(obj)
            Nk = ceil(obj.N * obj.CSCalculationMeshSizeCoef);
            [mu_temp,mu_weight_temp] = legzo_n1(Nk-2); 
            obj.mu_mesh_temp = [1 mu_temp -1]; 
            obj.mu_weight_temp = [0 mu_weight_temp 0];
        end
    end
    
    methods(Abstract = true, Access = protected)
        BaseCalculation(obj);
        AllocateResults(obj);
        AddExpansionCoef(obj,n,m,i,PlP);
        ExpandCrossSections(obj);
    end
    
end

