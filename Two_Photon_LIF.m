classdef Two_Photon_LIF %%< General_Operation_Point

    properties (Access=public)
        data_folder_path %Path to Folder containing .tif data for each timestep/axsial shift (string)
        data_folder_info %Struct containing name & path info about each timestep/axial shift data conatining folder (struct)
        data_info %Struct containing combined names of data containting folders and data itself (struct)
        
        data_quantity=0; %Quantity of .tif data for each timestep/axial shift (double)
        current_data %Changing Datastruct for each timestep/axial shift (struct)
        final_data %Final struct containing all necessary Parameters (struct) 
    end
    
    methods (Access=public)
        %Konstruktor sollte keine Incremente haben!!! "Structinput" muss in
        %Childclass selbst definiert und erweitert werden.
        function obj = Two_Photon_LIF()
            %obj@General_Operation_Point(Structinput,Varargin);
            obj.data_folder_path = uigetdir('Path to data containing folder');
            obj.data_folder_info = dir(obj.data_folder_path);
            %Initialisierung notwendig für mögliche Kombinierung der Structs 
            obj.data_info = dir([obj.data_folder_info(3).folder '\' obj.data_folder_info(3).name]);
            for k=4:numel(obj.data_folder_info)
                obj.data_info = [obj.data_info,dir([obj.data_folder_info(k).folder '\' obj.data_folder_info(k).name])];
            end
            for folder = 3:numel(obj.data_folder_info)
                for tif = 1:length(obj.data_info)
                    obj.current_data = Calculate_Data(obj.data_info,folder,tif);
                end
            end
        end
    end
    
    %events
        
    %end
    
    %enummeration
    
    %end
end

