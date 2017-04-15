classdef BaseLayerCalculation < BaseCalculation
    %BaseLayerCalculation Calculation for one layer
    %   Detailed explanation goes here
    
    properties (SetAccess = protected)
        Layer@Layer;                                    % Layer information and properties
    end
    
    properties (Dependent = true, SetAccess = protected)
        energy_mesh;                                    % Energy mesh
        Material;                                       % Material object in Layer shortcut
    end
    
    methods
        function obj = BaseLayerCalculation(Layer)
            
            assert(isa(Layer,'Layer'),'Layer must be of class Layer');
            
            obj = obj@BaseCalculation;
            
            obj.Layer = Layer;
            obj.ValidateMaterialE0;
            
            if ~obj.Layer.Material.IsManualDECS
              obj.Layer.Material.CalculateDECS;
            end
            
            addlistener(obj.Layer.Material,'EnergyChanged',@obj.ValidateMaterialE0);
            addlistener(obj.Layer.Material,'DECSChanged',@obj.ClearResults);
            addlistener(obj,'M','PostSet',@obj.ClearResults);
        end
        
        function val = get.energy_mesh(obj)
            val = obj.Layer.Material.DIIMFP_E;
        end
        
        function val = get.Material(obj)
            val = obj.Layer.Material;
        end
    end
    
    methods (Access = protected)
        function ValidateMaterialE0(obj,varargin)
            assert(isscalar(obj.Layer.Material.E0), 'Material E0 parameter must be scalar');
        end
    end
    
    methods (Abstract = true, Access = protected)
        CalculateCrossSection(obj,varargin);
    end
    
end

