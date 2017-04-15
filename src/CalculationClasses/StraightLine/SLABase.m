classdef SLABase < BaseLayerCalculation
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = protected)
        x_m_m                                   % Recalculated transmition DECS for all m values
        x_p_m                                   % Recalculated reflection DECS for all m values
    end
    
    methods
        function obj = SLABase(Layer)
            obj = obj@BaseLayerCalculation(Layer);

            obj.CalculateCrossSection;
        end
        
        function Calculate(obj)
            obj.ClearResults;
            obj.CalculateCrossSection;
            
            obj.BaseCalculation;
            
            obj.IsCalculated = true;
        end
    end
    methods(Access = protected)
        function CalculateCrossSection(obj,varargin)
%             [obj.x_p_m, obj.x_m_m] = X_mp(obj.Layer.Material.DECS,obj.Layer.Material.DECS_mu,obj.mu_mesh,obj.mu_mesh_weights,obj.M,1);
            x_l = obj.Layer.Material.Calculate_Leg_coefs(2000)*obj.norm;
            [obj.x_p_m, obj.x_m_m] = expandLegCoefsTo2ParamCrossSection(x_l, obj.mu_mesh,obj.M);
        end
    end
    
    methods(Abstract = true, Access = protected)
        BaseCalculation;
    end
    
end

