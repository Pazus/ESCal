classdef XPSAnode
    %XPSAnode A class to represent basic XPS anode properties
    
    properties (SetAccess = protected)
        MatName@char;                   % Material abbreviation
        PhotonEnergy@double scalar;     % Initial photon energy in eV
        polarization                    % photon polarization of X-Ray source
                                        % Polarization must be "unpolarized" or "circular" or "linear" or "none"
    end
    
    properties (Access = private, Constant = true);
        defaultPolarization = 'unpolarized'
        validPolarizations = {'unpolarized','circular','linear','none'}
    end
    
    properties (Access = protected, Constant = true)
        dataFileName = 'SourceData.mat';
    end
    
    methods
        function obj = XPSAnode(MatName, varargin)

            ps = inputParser;
            ps.FunctionName = 'XPSAnode';
            ps.addRequired('MatName',@(x) validateattributes(x,{'char'},{'nonempty'}));
            ps.addOptional('PhotonEnergy',[], @(x) validateattributes(x, {'numeric'},{'positive','scalar'}));
            ps.addOptional('Polarization',obj.defaultPolarization);

            ps.parse(MatName, varargin{:});
            
            validatestring(ps.Results.Polarization,obj.validPolarizations,'Polarization');

            obj.MatName = ps.Results.MatName;
            if isempty(ps.Results.PhotonEnergy)
               Data = obj.LoadDataFile(MatName);
               obj.PhotonEnergy = Data.PhotonEnergy;
            else
                obj.PhotonEnergy = ps.Results.PhotonEnergy;
            end;
            
            obj.polarization = ps.Results.Polarization;      
        end;
    end
    
    methods (Access = protected)
        function Data = LoadDataFile(obj, MatName)
            LoadData = load(obj.dataFileName);
            validatestring(MatName,fieldnames(LoadData.SourceData),'XPSAnode','MatName');
            Data = LoadData.SourceData.(MatName);
        end
    end;
    
end

