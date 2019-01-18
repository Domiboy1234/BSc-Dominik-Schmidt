 classdef General_Operation_Point
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = protected)
        folder % Place of OP data in the folder tree of the measurement
        current_path % Starting point for folder tree
        HeadFolder % Head of foldertree
        imagenames % Imagenames of the OP
        Pos_Nr % Number of sub OP e.g. Tiles or Points in time
        Copy_Template = 0;% Determines if Parameters for evaluation need to be updated
    end
    
    methods
        % Constructor  
        function obj = General_Operation_Point(strucinput,varargin)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            p = inputParser;          
            
              
            p.KeepUnmatched = true;
            addRequired(p,'strucinput',@isstruct);
            addParameter(p,'loaddata','No');
            defaulttemplate = struct;
            addParameter( p,'copytemplate',defaulttemplate);
            addParameter(p,'current_path','C:');
            parse(p,strucinput,varargin{:});
            
            strucinput = p.Results.strucinput;
            loaddata = p.Results.loaddata;
            copytemplate = p.Results.copytemplate;
            
            if strcmp(loaddata,'Yes')==1
                strucinput.current_path(:,1)  ={p.Results.current_path};
                obj = Import_Data(obj,strucinput);
                obj = check_path(obj);
            else
                if ~isempty(fieldnames(copytemplate))
                    obj = Import_Data(obj,copytemplate);
%                     obj.Copy_Template = 1;
                else
                    copytemplate.Copy_Template = 0;
                end
                obj.folder = strucinput.folder;
                obj.imagenames = strucinput.imagenames;
                obj.Pos_Nr = strucinput.Pos_Nr;
                obj.current_path = strucinput.current_path;
                obj.HeadFolder = strucinput.HeadFolder;
                obj = check_path(obj);
            end
        end
        
        function obj = check_path(obj)% Check if location of data is valid
            if isfolder([obj.current_path{1,1},'\',obj.folder{1,1}]) == 1
            else
                obj.current_path = uigetdir('C:\','Select current Data Path');
                % obj.current_path = {'H:'};
            end
        end
        
        % Import from structure
        function obj = Import_Data(obj,Data)
            propertylist = fieldnames(struct(obj));
            
            for i = 1 : length(propertylist)
                if isfield(Data,char(propertylist{i}))
                    obj.(char(propertylist{i})) = Data.(char(propertylist{i}));
                else
                    
                end
            end
        end
        
        % Export
        function Data = Export_Data(obj)
            Data = struct(obj);
            %             propertylist = properties(obj);
            %             for i = 1 : length(propertylist)
            %                 Data.(char(propertylist{i})) = obj.(char(propertylist{i}));
            %             end
        end
       function [Switch,argout] = Filter(obj,varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            addRequired(p,'obj');
            addParameter(p,'Name_Search',{'default'});
            addParameter(p,'Value_Search',{'default'});
            parse(p,obj,varargin{:});
            obj = p.Results.obj;
            Name_Search = p.Results.Name_Search;
            Value_Search = p.Results.Value_Search;
            
            if strcmp(Name_Search{1},'default')
                output_ok = 0;
                found_error = 0;
                while output_ok == 0
                    if found_error == 1
                        uiwait(msgbox('Number of Names does not match Number of Values', 'Error','error'));
                    end
                    prompt = {'Enter Names for Operation Condition seperated by ";"','Enter Values for Operation Condition seperated by ";"'};
                    dlg_title = 'Input Filter';
                    f1 = inputdlg(prompt,dlg_title,1,{strjoin(obj.Docu_Names,';'),strjoin(obj.Docu_Values,';')});
                    Name_Search = strsplit(f1{1},';');
                    length_Names = length(obj.Docu_Names);
                    
                    Value_Search = strsplit(f1{2},';');
                    length_Values = length(obj.Docu_Values);
                    if  length_Names==length_Values
                        output_ok = 1;
                    else
                        found_error = 1;
                    end
                end
            end
            
            f = zeros(length(Name_Search),1);
            for k = 1 : length(f)
                %         f0 = find(contains(Condition_Names,searchstr_Name{k}));
                f0 = strfind(obj.Docu_Names,Name_Search{k});
                f0 = find(not(cellfun('isempty', f0)));
                if strcmp(obj.Docu_Values{f0},Value_Search{k})==1
                    f(k)=1;
                end
            end
            if sum(f) == length(f)
                Switch = 1;
            else
                Switch = 0;   
            end
            argout{1} = 'Name_Search';
            argout{2} = Name_Search;
            argout{3} = 'Value_Search';
            argout{4} = Value_Search;
        end 
    end
    
%     methods (Abstract)
%         Evaluation(obj)
%     end
end

