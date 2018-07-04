classdef MaterialPropertiesDatabase

    methods (Static)
        function data=getData(materialName)
            persistent databaseCache;
            
            if isempty(databaseCache) || ~isfield(databaseCache,materialName)
                MaterialDatabase = load('MaterialData.mat','MaterialData');
                validatestring(materialName,fieldnames(MaterialDatabase.MaterialData),'Material','materialName');
                databaseCache.(materialName) = MaterialDatabase.MaterialData.(materialName);
            end
            
            data=databaseCache.(materialName);
        end
        
        function list = getMaterialList()

                MaterialDatabase = load('MaterialData.mat','MaterialData');
                list = fields(MaterialDatabase.MaterialData);

        end
    end
    
end

