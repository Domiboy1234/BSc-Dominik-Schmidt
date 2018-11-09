classdef Highspeed_Operation_Point < General_Operation_Point
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Hidden)
        frequency % frequency of data aquisition
        Numberofpics % Number of pictures for a single Experiment
        BackgroundSequence % Index for pictures for Backgroundcorrection in a single shot
        NumberofBackgroundpics % Number of pictures for Backgroundcorrection in a single shot
        scale % scale in mm/px
        Date
        Clocktime
        time
    end
    
    methods
        
        % Constructor
        function obj = Highspeed_Operation_Point(strucinput,varargin)% Create OP with class HS_XRay_MultiExposure and define picture parameters for evaluation
            obj@General_Operation_Point(strucinput,varargin{:});
             p = inputParser;
             p.KeepUnmatched = true;
            addRequired(p,'strucinput',@isstruct);
            addParameter(p,'loaddata','No');
            parse(p,strucinput,varargin{:});
%              strucinput = p.Results.strucinput;
            loaddata = p.Results.loaddata;
            
            if ~strcmp(loaddata,'Yes')==1
                obj = get_cih_data(obj);
                if obj.Copy_Template == 0
                    obj = def_scale(obj);
                    obj = def_HSparameters(obj);
                end
            end
        end
        % HS Parameters
        function obj = get_cih_data(obj)
            obj = check_path(obj);
            cih_file = dir([obj.current_path{1,1},'\',obj.folder{1,1},'\*.cih']);
            output = cih_import_func([cih_file.folder,'\',cih_file.name]);
            obj.Date = char(output{2,3});
            obj.Clocktime = char(output{3,3});
            obj.frequency = str2double(output{16,4})/1000;
        end
        function obj = def_scale(obj)% define scale
            output_ok = 0;
            found_error = 0;
            while output_ok == 0
                if found_error == 1
                    uiwait(msgbox('No scale Parameter available from workspace, please try again!', 'Error','error'));
                end
                [file_scale, path_scale] = uigetfile({'*.tif;*.tiff;*.jpg;*.bmp','All Image Files'},'Select original image for scale definition, cancel if you want to load the value from workspace');
                if ~ischar(file_scale)
                    nopic = 1;
                else
                    pic = imread([path_scale '\' file_scale]);
                    pic=double(pic);
                    pic=pic./max(max(pic));
                    nopic = 0;
                end
                if nopic == 1 && isempty(obj.scale)==0
                    obj.scale = obj.scale;
                    output_ok = 1;
                elseif nopic == 0
                    figure('NumberTitle','off','Name','Definition of scale using two points.');
                    imshow(pic);
                    title('Press enter to select two points after zooming!');
                    zoom on;
                    waitfor(gcf,'CurrentCharacter',char(13));
                    zoom reset;
                    zoom off;
                    [x_y_scale(1,:)] = ginput(1); hold on
                    plot(x_y_scale(1,1), x_y_scale(1,2), 'rx', 'MarkerSize', 10);
                    [x_y_scale(2,:)] = ginput(1); hold on
                    plot(x_y_scale(2,1), x_y_scale(2,2), 'rx', 'MarkerSize', 10);
                    pause(1);
                    close(gcf);
                    pixels = ((x_y_scale(1,1)-x_y_scale(2,1))^2 + (x_y_scale(1,2)-x_y_scale(2,2))^2)^0.5;
                    
                    prompt = {'Enter distance in mm:'};
                    dlg_title = 'Input';
                    distance = {inputdlg(prompt, dlg_title)};
                    obj.scale = str2double(distance{1,1}) / pixels;
                    output_ok = 1;
                else
                    found_error = 1;
                end
            end
        end
        function obj = def_HSparameters(obj)% define parameters for highspeed evaluation
            output_ok = 0;
            found_error = 0;
            while output_ok == 0
                if found_error == 1
                    uiwait(msgbox('Incorrect Input, please try again!', 'Error','error'));
                end
                prompt = {'Enter frequency in kHz:','Enter Number of Pictures during a single event:','Enter Sequence of Background pictures during a single event:'};
                f= inputdlg(prompt,'Input',1,{num2str(obj.frequency),num2str(obj.Numberofpics),num2str(obj.NumberofBackgroundpics)});
                if isempty(f{1,1})==1 || isempty(f{2,1})==1 || isempty(f{3,1})==1
                else
                    obj.frequency = str2double(f{1});
                    obj.Numberofpics = str2double(f{2});
                    sequence = strsplit(f{3},':');
                    obj.BackgroundSequence = str2double(sequence(1)):1:str2double(sequence(length(sequence)));
                    obj.NumberofBackgroundpics = length(obj.BackgroundSequence);
                    if isnan(obj.frequency) == 0 && isnan(obj.Numberofpics) == 0 && isnan(obj.NumberofBackgroundpics) == 0
                        output_ok = 1;
                    end
                end
                found_error = 1;
            end
        end
        function time = get.time(obj)
            time = 0:10^6 / (1000*obj.frequency):(10^6 / (1000*obj.frequency)*obj.Numberofpics-10^6 / (1000*obj.frequency));
        end

    end
    methods (Abstract)
        Evaluation(obj)

    end
end

