classdef (Abstract)  CoreMaterial < matlab.mixin.Copyable
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetObservable = true, AbortSet = true)
        E0;
        useTransportApproximation = false;
    end
    
    properties (SetAccess = protected)
        M;
        Z;
        Density;

        lambda;
        
        l_in;
        l_el;
        l_tr;
        l_tot;
        
        DIIMFP;
        DIIMFP_E;
        DIIMFP_Shielding_koef;
        IsManualDIIMFP@logical scalar = false;
        
        DECS; % Differential elastic cross section
        DECS_mu; % Mesh on which the differential elastic cross section is calculated
        DECS_l; % Legendre polynomial expansion of DECS
        
        IsManualDECS@logical scalar = false;
        IsManualIMFP@logical scalar = false;
        IsManualEMFP@logical scalar = false;
        IsManualTRMFP@logical scalar = false;
    end
    
    properties (Access = private, Constant = true)
        min_DECS_mu_mesh_size = 50;
    end
    
    properties (Dependent = true, SetAccess = protected)
        sigma_el
        sigma_in
        sigma_tr
    end
    
    events
        EnergyChanged;
        DIIMFPChanged;
        DECSChanged;
    end
    
    methods
        function RecalculateAll(obj, varargin)   
            obj.CalculateAllMFP;
            if ~isempty(obj.DECS) 
                if obj.IsManualDECS
                    disp('DECS is manually set and was not recalculated.');
                else
                    obj.CalculateDECS(obj.DECS_mu);
                end;
                 
            end
            if ~isempty(obj.DIIMFP)
                if obj.IsManualDIIMFP
                    disp('DIIMFP is manually set and was not recalculated.');
                else
                    obj.CalculateDIIMFP(obj.DIIMFP_E,obj.DIIMFP_Shielding_koef);
                end;
            end;
            notify(obj,'EnergyChanged');
        end
        
        function CalculateAllMFP(obj)
            if obj.IsManualIMFP
                disp('IMFP is manually set and was not recalculated.');
            else
                obj.CalculateIMFP;
            end;
            if obj.IsManualEMFP
                disp('EMFP is manually set and was not recalculated.');
            else
                obj.CalculateEMFP;
            end;   
            if obj.IsManualTRMFP
                disp('TRMFP is manually set and was not recalculated.');
            else
                obj.CalculateTRMFP;
            end;
            
            obj.l_tot = obj.l_in.*obj.l_el./(obj.l_in + obj.l_el);
            obj.lambda = obj.l_tot./obj.l_el;
        end
        
        function SetManualDECS(obj,mesh_mu,DECS)
            
            validateattributes(mesh_mu,{'numeric'},{'nonempty','vector','column'});
            validateattributes(DECS,{'numeric'},{'nonempty','nonnegative','size',[numel(mesh_mu), numel(obj.E0)]});
            
            assert(numel(mesh_mu) >= obj.min_DECS_mu_mesh_size, ...
                'Number of elements of mu_mesh has to be more or equal to %d', ...
                obj.min_DECS_mu_mesh_size);
            
            
            if max(mesh_mu)<=1 && min(mesh_mu)>=-1
                mesh_mu = acos(mesh_mu);
            elseif max(mesh_mu)>pi && min(mesh_mu)>=0
                mesh_mu = mesh_mu/180*pi;
            end
            
            obj.DECS = DECS;
            obj.DECS_mu = reshape(mesh_mu,[],1);
            
            obj.IsManualDECS = true;
            
            notify(obj,'DECSChanged');
        end
        
        function SetManualDIIMFP(obj, mesh_E, DIIMFP)
            
            
            validateattributes(mesh_E,{'numeric'},{'nonempty','vector', 'column', 'nonnegative'}, 'SetManualDIIMFP', 'mesh_E');
%             validateattributes(DIIMFP,{'numeric'},{'nonempty','nonnegative','size',[numel(mesh_E), numel(obj.E0)]}, 'SetManualDIIMFP', 'DIIMFP');
            
            if mesh_E(2)-mesh_E(1)<0
                error('Check the energy mesh for DIIMFP: energy values must rise from min to max!');
            end
            
            obj.DIIMFP = DIIMFP;
            obj.DIIMFP_E = reshape(mesh_E,[],1);
            
            obj.IsManualDIIMFP = true;
            
            notify(obj,'DIIMFPChanged');
        end
        
        function SetManualIMFP(obj, IMFP)
            validateattributes(IMFP,{'numeric'},{'nonempty','column','nonnegative','size',[numel(obj.E0) 1]});
            
            obj.l_in = IMFP;
            
            obj.IsManualIMFP = true;
        end
        
        function SetManualEMFP(obj, EMFP)
            persistent ps;
            if isempty(ps)
                ps = inputParser;
                ps.FunctionName = 'SetManualEMFP';
                ps.addRequired('EMFP',@(x) validateattributes(x,{'numeric'},{'nonempty','column','nonnegative','size',[numel(obj.E0) 1]}));
            end
            
            ps.parse(EMFP);
            
            obj.l_el = EMFP;
            
            obj.IsManualEMFP = true;
        end
        
        function SetManualTRMFP(obj, TRMFP)
            persistent ps;
            if isempty(ps)
                ps = inputParser;
                ps.FunctionName = 'SetManualTRMFP';
                ps.addRequired('TRMFP',@(x) validateattributes(x,{'numeric'},{'nonempty','column','nonnegative','size',[numel(obj.E0) 1]}));
            end
            
            ps.parse(TRMFP);
            
            obj.l_el = TRMFP;
            
            obj.IsManualTRMFP = true;
        end
        
        function val = get.sigma_el(obj)
            if isempty(obj.l_el)
                val = [];
            else
                val = 1./(obj.Density*obj.l_el);
            end
        end
        function val = get.sigma_in(obj)
            if isempty(obj.l_el)
                val = [];
            else
                val = 1./(obj.Density*obj.l_in);
            end
        end
        function val = get.sigma_tr(obj)
            if isempty(obj.l_el)
                val = [];
            else
                val = 1./(obj.Density*obj.l_tr);
            end
        end
        
        function plotDECS(obj)
            if isempty(obj.DECS); disp('DECS cross section was not yet calculated'); return; end
            semilogy(obj.DECS_mu,obj.DECS);
            xlabel('\mu');
            ylabel('x_{el}(\mu)');
        end
        
        function plotDIIMFP(obj)
            if isempty(obj.DIIMFP); disp('DIIMFP cross section was not yet calculated'); return; end
            
            plot(obj.DIIMFP_E,obj.DIIMFP);
            xlabel('E')
            ylabel('x_{in}(E)');
            if numel(obj.E0) < 5
                legend_str = cell(1,numel(obj.E0));
                for i=1:numel(obj.E0)
                    legend_str{i} = ['E0 = ' num2str(obj.E0(i))];
                end
                legend(legend_str);
            end
        end
        
    end
    
    methods (Abstract)
        [DECS_temp] = CalculateDECS(obj, mu);
        
        CalculateDIIMFP(obj, varargin);
        
        Res =  isequal(mat1, mat2);
    end
    
    methods (Abstract, Access = protected)
        l_in = CalculateIMFP(obj);
        l_el = CalculateEMFP(obj);
        l_tr = CalculateTRMFP(obj);
    end
    
end

