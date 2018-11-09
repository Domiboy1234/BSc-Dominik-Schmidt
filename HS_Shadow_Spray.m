classdef HS_Shadow_Spray
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    properties (Hidden)
        % Parameters for Evaluation
        folder % Place of OP data in the folder tree of the measurement
        current_path % Starting point for folder tree
        imagenames % Imagenames of the OP
        Pos_Nr % Number of sub OP e.g. Tiles or Points in time
        scale % scale in mm/px
        position_nozzle % y and x coordinate of the nozzle tip
        bw_mask_spray % mask for boundary definition bw picture
        poly_mask_spray % mask for boundary definition
        frequency % frequency of data aquisition
        Numberofpics % Number of pictures for a single Experiment
        NumberofBackgroundpics % Number of pictures for Backgroundcorrection in a single shot
        threshold % Threshold for Background correction
        bwarea %
        bw_mask_nozzle % mask for lowest pixel value bw picture
        poly_mask_nozzle % mask for lowest pixel value
        Docu_Names % Documentation Variables
        Docu_Values % Documentation Values
        Copy_Template = 0;% Determines if Parameters for evaluation need to be updated
        sprayXmask_crossing % Determines whether the spray has reached the edge of the detectable area
        Date
        Clocktime
    end
    properties (Dependent,Access = private)
        time
    end
    properties
        % Parameters from Evaluation
        axial_penetration
        total_penetration
        plumeangle_2planes
        plumeangle_maxspraywidth
        plumeangle_circle
        spraydirection
        sprayarea
    end
    
    methods
        % Constructor
        function obj = HS_Shadow_Spray(strucinput,varargin)% Create OP with class HS_Shadow_Spray and define picture parameters for evaluation
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            %             strucinput = varargin{1};
            %             if nargin > 1
            %                 obj = varargin{2};
            %             end
            p = inputParser;
            addRequired(p,'strucinput',@isstruct);
            addParameter(p,'loaddata','No');
            defaulttemplate = struct;
            addParameter(p,'copytemplate',defaulttemplate);
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
                else
                    copytemplate.Copy_Template=0;
                end
                obj.folder = strucinput.folder;
                obj.imagenames = strucinput.imagenames;
                obj.Pos_Nr = strucinput.Pos_Nr;
                obj.current_path = strucinput.current_path;
                obj = check_path(obj);
                if ~isempty(copytemplate) && copytemplate.Copy_Template==1
                else % Acquisition of Parameters for Evaluation
                    obj = def_scale(obj);
                    obj = def_nozzle(obj);
                    obj = def_spraymask(obj);
                    obj = def_HSparameters(obj);
                    obj = def_BWparameters(obj);
                    obj = def_nozzlemask(obj);
                    obj = get_cih_data(obj);
                end
                obj = def_OPdocu(obj);
                obj.Copy_Template = 0;
            end
        end
        % Parameter input, happens in the Constructor
        function obj = check_path(obj)% Check if location of data is valid
            if isfolder([obj.current_path{1,1},'\',obj.folder{1,1}]) == 1
            else
                cp = uigetdir('C:\','Select current Data Path');
                obj.current_path = {cp(1:length(cp)-1)};
%                 obj.current_path = {'F:'};
            end
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
        function obj = def_nozzle(obj) % define nozzle position
            obj = check_path(obj);
            picture = imread([obj.current_path{1,1},'\',obj.folder{1,1},'\',obj.imagenames{1,2}]);
            picture = double(picture)/max(max(double(picture)));
            loadpara = 'No';
            imshow(picture);
            title('Press enter to select a point for nozzle definition after zooming!');
            if isempty(obj.position_nozzle)==0
                hold on
                plot(obj.position_nozzle(1,2), obj.position_nozzle(1,1), 'go', 'MarkerSize', 10);
                loadpara = questdlg('Use previous coordinates for nozzle?','Input', 'Yes', 'No','Yes');
            end
            
            if strcmp(loadpara,'No')==1
                zoom on;
                waitfor(gcf,'CurrentCharacter',char(13));
                zoom reset;
                zoom off;
                [obj.position_nozzle(1,2), obj.position_nozzle(1,1)] = ginput(1); hold on
            elseif strcmp(loadpara,'Yes')==1
                obj.position_nozzle = obj.position_nozzle;
            end
            plot(obj.position_nozzle(1,2), obj.position_nozzle(1,1), 'rx', 'MarkerSize', 10);
            pause(1);
            close(gcf);
        end
        function obj = def_spraymask(obj)% define mask for maximum spray propagation
            obj = check_path(obj);
            picture = imread([obj.current_path{1,1},'\',obj.folder{1,1},'\',obj.imagenames{1,2}]);
            picture = double(picture)/max(max(double(picture)));
            figure('NumberTitle','off','Name','Select mask for spray');
            imshow(picture);
            if isempty(obj.poly_mask_spray) == 0
                h = impoly(gca,obj.poly_mask_spray);
            else
                h = impoly(gca);
            end
            obj.poly_mask_spray = wait(h);
            close(gcf);
            obj.bw_mask_spray = poly2mask(obj.poly_mask_spray(:,1), obj.poly_mask_spray(:,2), size(picture,1), size(picture,2));
        end
        function obj = def_HSparameters(obj)% define parameters for highspeed evaluation
            output_ok = 0;
            found_error = 0;
            while output_ok == 0
                if found_error == 1
                    uiwait(msgbox('Incorrect Input, please try again!', 'Error','error'));
                end
                prompt = {'Enter frequency in kHz:','Enter Number of Pictures during a single event:','Enter Number of Background pictures during a single event:'};
                f= inputdlg(prompt,'Input',1,{num2str(obj.frequency),num2str(obj.Numberofpics),num2str(obj.NumberofBackgroundpics)});
                if isempty(f{1,1})==1 || isempty(f{2,1})==1 || isempty(f{3,1})==1
                else
                    obj.frequency = str2double(f{1});
                    obj.Numberofpics = str2double(f{2});
                    obj.NumberofBackgroundpics = str2double(f{3});
                    if isnan(obj.frequency) == 0 && isnan(obj.Numberofpics) == 0 && isnan(obj.NumberofBackgroundpics) == 0
                        output_ok = 1;
                    end
                end
                found_error = 1;
            end
        end
        function obj = def_BWparameters(obj)% define parameters for binarization
            output_ok = 0;
            found_error = 0;
            while output_ok == 0
                if found_error == 1
                    uiwait(msgbox('Incorrect Input, please try again!', 'Error','error'));
                end
                prompt = {'Enter Thresholdvalue between 0 and 1 or "auto" for automatic Threshold:','Enter Value for BWarea:'};
                f= inputdlg(prompt,'Input',1,{num2str(obj.threshold),num2str(obj.bwarea)});
                if isempty(f{1,1})==1 ||  isempty(f{2,1})==1
                else
                    obj.threshold=str2double(f{1});
                    obj.bwarea=str2double(f{2});
                    if isnan(obj.threshold) == 0 && isnan(obj.bwarea) == 0
                        output_ok = 1;
                    end
                end
                found_error = 1;
            end
        end
        function obj = def_nozzlemask(obj)% define mask for image correction, darkest area is the nozzle tip
            obj = check_path(obj);
            picture = imread([obj.current_path{1,1},'\',obj.folder{1,1},'\',obj.imagenames{1,2}]);
            picture = double(picture)/max(max(double(picture)));
            figure('NumberTitle','off','Name','Select mask on nozzle for brightness adjustment');
            imshow(picture);
            zoom on;
            waitfor(gcf,'CurrentCharacter',char(13));
            zoom reset;
            zoom off;
            if isempty(obj.poly_mask_nozzle) == 0
                h = imrect(gca,obj.poly_mask_nozzle);
            else
                h = imrect(gca);
            end
            obj.poly_mask_nozzle = wait(h);
            close(gcf);
            x_mask = [obj.poly_mask_nozzle(1), obj.poly_mask_nozzle(1)+obj.poly_mask_nozzle(3), obj.poly_mask_nozzle(1)+obj.poly_mask_nozzle(3), obj.poly_mask_nozzle(1)];
            y_mask = [obj.poly_mask_nozzle(2), obj.poly_mask_nozzle(2), obj.poly_mask_nozzle(2)+obj.poly_mask_nozzle(4), obj.poly_mask_nozzle(2)+obj.poly_mask_nozzle(4)];
            obj.bw_mask_nozzle = poly2mask(x_mask, y_mask, size(picture,1), size(picture,2));
        end
        function obj = def_OPdocu(obj)% documentation of OP
            found_error = 0;
            output_ok = 0;
            while output_ok == 0
                if found_error == 1
                    uiwait(msgbox('Number of Names does not match Number of Values', 'Error','error'));
                end
                prompt = {'Enter Names for Operation Condition seperated by ";"','Enter Values for Operation Condition seperated by ";"'};
                dlg_title = 'Input';
                if isempty(obj.Docu_Names) == 0
                    f1 = inputdlg(prompt,dlg_title,1,{strjoin(obj.Docu_Names,';'),strjoin(obj.Docu_Values,';')});
                else
                    f1 = inputdlg(prompt,dlg_title);
                end
                if ~isempty(f1)
                    obj.Docu_Names = strsplit(f1{1},';');
                    length_Names = length(obj.Docu_Names);
                    
                    obj.Docu_Values = strsplit(f1{2},';');
                    length_Values = length(obj.Docu_Values);
                    if  length_Names==length_Values
                        output_ok = 1;
                    end
                end
                found_error = 1;
            end
        end
        function obj = get_cih_data(obj)
            obj = check_path(obj);
            cih_file = dir([obj.current_path{1,1},'\',obj.folder{1,1},'\*.cih']);
            output = cih_import_func([cih_file.folder,'\',cih_file.name]);
            obj.Date = char(output{2,3});
            obj.Clocktime = char(output{3,3});
        end
        
        % Data Evaluation
        
        function obj = Evaluation(obj)
            obj = check_path(obj);
            obj = initialize_output(obj);
            obj = get_cih_data(obj);
            for i = 1:size(obj.imagenames,2)
                seriesnumb = ceil(i/obj.Numberofpics);
                imagenumb = i-(seriesnumb-1)*obj.Numberofpics;
                if imagenumb == 1
                    back_mean_normed = eval_background(obj,i);
                end
                [picture_normed,picture_full_adjust_imcom] = eval_imageadjust(obj,i,back_mean_normed);
                [picture_edge, picture_bwarea] = eval_binarize(obj,picture_full_adjust_imcom);
                obj = eval_sprayXmask_crossing(obj,i,picture_edge);
                obj = eval_axial_penetration(obj,i,picture_edge);
                obj = eval_total_penetration(obj,i,picture_edge);
                obj = eval_plumeangle_2planes(obj,i,picture_edge);
                obj = eval_plumeangle_maxspraywidth(obj,i,picture_edge,picture_bwarea,picture_full_adjust_imcom);
                obj = eval_spraydirection(obj,i,picture_edge,picture_bwarea,picture_full_adjust_imcom);
                obj = eval_plumeangle_circle(obj,i,picture_edge);
                obj = eval_sprayarea(obj,i,picture_edge,picture_bwarea);
%                 Export_ImageEdge(obj,i,picture_edge,picture_normed);
            end
        end
        function obj = initialize_output(obj)
            imagenumbmax = obj.Numberofpics;
            seriesnumbmax = floor(size(obj.imagenames,2)/obj.Numberofpics);
            %             obj.time = 0:10^6 / (1000*obj.frequency):(10^6 / (1000*obj.frequency)*obj.Numberofpics-10^6 / (1000*obj.frequency));
            obj.axial_penetration = NaN(imagenumbmax,seriesnumbmax);
            obj.total_penetration = NaN(imagenumbmax,seriesnumbmax);
            obj.plumeangle_2planes = NaN(imagenumbmax,seriesnumbmax);
            obj.plumeangle_maxspraywidth = NaN(imagenumbmax,seriesnumbmax);
            obj.plumeangle_circle = NaN(imagenumbmax,seriesnumbmax);
            obj.spraydirection = NaN(imagenumbmax,seriesnumbmax);
            obj.sprayarea = NaN(imagenumbmax,seriesnumbmax);
            obj.sprayXmask_crossing = ones(imagenumbmax,seriesnumbmax);
        end
        function back_mean_normed = eval_background(obj,imagenr)
            obj = check_path(obj);
            seriesnumb = ceil(imagenr/obj.Numberofpics);
            imagenumb = imagenr-(seriesnumb-1)*obj.Numberofpics;
            back_mean_normed = zeros(size(obj.bw_mask_spray,1),size(obj.bw_mask_spray,2),obj.NumberofBackgroundpics);
            j = 1;
            for i = 1 + obj.Numberofpics * (seriesnumb-1) : obj.NumberofBackgroundpics + obj.Numberofpics * (seriesnumb-1)
                background = imread([obj.current_path{1,1},'\',obj.folder{1,1},'\',obj.imagenames{1,i}]);
                back_mean_normed(:,:,j) = double(background) / max(max(double(background).*obj.bw_mask_spray));
                j = j+1;
            end
            back_mean_normed(back_mean_normed>1)=1;
            back_mean_normed = mean(back_mean_normed,3);
        end
        function [picture_normed,picture_full_adjust_imcom] = eval_imageadjust(obj,imagenr,back_mean_normed)
            obj = check_path(obj);
            picture = imread([obj.current_path{1,1},'\',obj.folder{1,1},'\',obj.imagenames{1,imagenr}]);
            picture_normed = double(picture)/max(max(double(picture).*obj.bw_mask_spray));
            picture_normed(picture_normed>1)=1;
            seriesnumb = ceil(imagenr/obj.Numberofpics);
            imagenumb = imagenr-(seriesnumb-1)*obj.Numberofpics;
            black_nozzle = mean(picture_normed(obj.bw_mask_nozzle));
            picture_adjust = imadjust(picture_normed, [black_nozzle;1], [0;1]);
            black_nozzle = mean(back_mean_normed(obj.bw_mask_nozzle));
            back_mean_adjust = imadjust(back_mean_normed, [black_nozzle;1], [0;1]);
            picture_back_corrected = back_mean_adjust - picture_adjust;
            bw_back_corrected = imbinarize(picture_back_corrected, obj.threshold);
            bw_back_corrected = bwareaopen(bw_back_corrected, obj.bwarea);
            picture_adjust_masked = picture_adjust;
            picture_adjust_masked(~bw_back_corrected) = 0;
            masked_brightest_point= max(max(picture_adjust_masked));
            if masked_brightest_point == 0
                picture_full_adjust_imcom = zeros(size(picture_adjust));
            else
                picture_full_adjust = imadjust(picture_adjust, [0; masked_brightest_point], [0;1]);
                picture_full_adjust(~bw_back_corrected) = 1;
                picture_full_adjust_imcom = imcomplement(picture_full_adjust);
            end
        end
        function [picture_edge,picture_bwarea] = eval_binarize(obj,picture_full_adjust_imcom)
            picture_thresh = imbinarize(picture_full_adjust_imcom.*obj.bw_mask_spray, obj.threshold);
            picture_bwarea = bwareaopen(picture_thresh, obj.bwarea);
            picture_bwarea = imcomplement(picture_bwarea);
            picture_bwarea = bwareaopen(picture_bwarea,obj.bwarea);
            picture_bwarea = imcomplement(picture_bwarea);
            picture_edge = edge(picture_bwarea, 'canny');
        end
        
        function obj = eval_sprayXmask_crossing(obj,imagenr,picture_edge)
            center_nozzle_vec = [size(obj.bw_mask_spray,1)/2 - obj.position_nozzle(1,1),size(obj.bw_mask_spray,2)/2 - obj.position_nozzle(1,2)];
            seriesnumb = ceil(imagenr/obj.Numberofpics);
            imagenumb = imagenr-(seriesnumb-1)*obj.Numberofpics;
            circle_radius = 5;
            sprayXmask_crossing_picture = picture_edge .* edge(obj.bw_mask_spray, 'canny');
            [row, col] = find(sprayXmask_crossing_picture); %row = y, col = x
            distance = round((( obj.position_nozzle(1,1)-row).^2+( obj.position_nozzle(1,2)-col).^2).^0.5);
            d = find((round(circle_radius/obj.scale) - 1) >= distance);
            sprayXmask_crossing_picture(row(d),col(d))=0;
            if sum(sum(sprayXmask_crossing_picture))>0
                obj.sprayXmask_crossing(imagenumb,seriesnumb) = NaN;
            end
        end
        function obj = eval_axial_penetration(obj,imagenr,picture_edge)
            center_nozzle_vec = [size(obj.bw_mask_spray,1)/2 - obj.position_nozzle(1,1),size(obj.bw_mask_spray,2)/2 - obj.position_nozzle(1,2)];
            seriesnumb = ceil(imagenr/obj.Numberofpics);
            imagenumb = imagenr-(seriesnumb-1)*obj.Numberofpics;
            [row, col] = find(picture_edge); %row = y, col = x
            edge_index = [row, col];
            f_coincidential = find(abs(center_nozzle_vec)==max(abs(center_nozzle_vec)));
            if isempty(row)
                obj.axial_penetration(imagenumb,seriesnumb) = NaN;
            else
                obj.axial_penetration(imagenumb,seriesnumb) = max(abs(edge_index(f_coincidential)-obj.position_nozzle(f_coincidential))) * obj.scale;
                seriesnumbmax = floor(size(obj.imagenames,2)/obj.Numberofpics);
                if sum(isnan(obj.sprayXmask_crossing(imagenumb,:)))> floor(0.2*seriesnumbmax)
                    obj.axial_penetration(imagenumb,:) = NaN;
                else
                    sprayXvector_crossing = mean(obj.sprayXmask_crossing(imagenumb,seriesnumb));
                    obj.axial_penetration(imagenumb,seriesnumb) = obj.axial_penetration(imagenumb,seriesnumb)*sprayXvector_crossing;
                end
            end
        end
        function obj = eval_total_penetration(obj,imagenr,picture_edge)
            %             center_nozzle_vec = [size(obj.bw_mask_spray,1)/2 - obj.position_nozzle(1,1),size(obj.bw_mask_spray,2)/2 - obj.position_nozzle(1,2)];
            seriesnumb = ceil(imagenr/obj.Numberofpics);
            imagenumb = imagenr-(seriesnumb-1)*obj.Numberofpics;
            [row, col] = find(picture_edge); %row = y, col = x
            %             edge_index = [row, col];
            %             f_coincidential = find(abs(center_nozzle_vec)==max(abs(center_nozzle_vec)));
            if isempty(row)
                obj.total_penetration(imagenumb,seriesnumb) = NaN;
            else
                obj.total_penetration(imagenumb,seriesnumb) = max(((obj.position_nozzle(1,1)-row).^2 + (obj.position_nozzle(1,2)-col).^2).^0.5) * obj.scale;
                seriesnumbmax = floor(size(obj.imagenames,2)/obj.Numberofpics);
                if seriesnumb == seriesnumbmax
                    test=1;
                end
                if sum(isnan(obj.sprayXmask_crossing(imagenumb,:)))> floor(0.2*seriesnumbmax)
                    obj.total_penetration(imagenumb,:) = NaN;
                else
                    sprayXvector_crossing = mean(obj.sprayXmask_crossing(imagenumb,seriesnumb));
                    obj.total_penetration(imagenumb,seriesnumb) = obj.total_penetration(imagenumb,seriesnumb)*sprayXvector_crossing;
                end

            end
        end
        function obj = eval_plumeangle_2planes(obj,imagenr,picture_edge)
            center_nozzle_vec = [size(obj.bw_mask_spray,1)/2 - obj.position_nozzle(1,1),size(obj.bw_mask_spray,2)/2 - obj.position_nozzle(1,2)];
            seriesnumb = ceil(imagenr/obj.Numberofpics);
            imagenumb = imagenr-(seriesnumb-1)*obj.Numberofpics;
            [row, col] = find(picture_edge); %row = y, col = x
            edge_index = [row, col];
            f_coincidential = find(abs(center_nozzle_vec)==max(abs(center_nozzle_vec)));
            if isempty(row)
                obj.plumeangle_2planes(imagenumb,seriesnumb) = NaN;
            else
                distance_1 = round(obj.position_nozzle(f_coincidential) +sign(center_nozzle_vec(f_coincidential)) * 5/obj.scale);
                distance_2 = round(obj.position_nozzle(f_coincidential) +sign(center_nozzle_vec(f_coincidential)) * 15/obj.scale);
                
                switch f_coincidential
                    case 2
                        y1_1 = max(find(picture_edge(:,distance_1)));
                        y2_1 = min(find(picture_edge(:,distance_1)));
                        y1_2 = max(find(picture_edge(:,distance_2)));
                        y2_2 = min(find(picture_edge(:,distance_2)));
                    case 1
                        y1_1 = max(find(picture_edge(distance_1,:)));
                        y2_1 = min(find(picture_edge(distance_1,:)));
                        y1_2 = max(find(picture_edge(distance_2,:)));
                        y2_2 = min(find(picture_edge(distance_2,:)));
                end
                
                
                vect1_ = [(distance_2-distance_1) (y1_2-y1_1)];
                vect2_ = [(distance_2-distance_1) (y2_2-y2_1)];
                obj.plumeangle_2planes(imagenumb,seriesnumb) = acos(dot(vect1_, vect2_)/sqrt(sum(vect1_.^2)*sum(vect2_.^2)))*180/pi;
            end
        end
        function obj = eval_plumeangle_maxspraywidth(obj,imagenr,picture_edge,picture_bwarea,picture_full_adjust_imcom)
            center_nozzle_vec = [size(obj.bw_mask_spray,1)/2 - obj.position_nozzle(1,1),size(obj.bw_mask_spray,2)/2 - obj.position_nozzle(1,2)];
            seriesnumb = ceil(imagenr/obj.Numberofpics);
            imagenumb = imagenr-(seriesnumb-1)*obj.Numberofpics;
            [row, col] = find(picture_edge); %row = y, col = x
            edge_index = [row, col];
            f_coincidential = find(abs(center_nozzle_vec)==max(abs(center_nozzle_vec)));
            if f_coincidential == 1
                f_perpendicular = 2;
            elseif f_coincidential == 2
                f_perpendicular = 1;
            end
            if isempty(row)
                plume_angle_max_spray_width1 = NaN;
                plume_angle_max_spray_width2 = NaN;
                plume_angle_max_spray_width3 = NaN;
                plume_angle_max_spray_width4 = NaN;
                plume_angle_max_spray_width5 = NaN;
                plume_angle_max_spray_width6 = NaN;
            else
                picture_bw_label = bwlabel(picture_bwarea);
                S = (struct2vektorWC(regionprops(picture_bw_label, picture_full_adjust_imcom,'WeightedCentroid')));
                
                Sc(1,2) = mean(S(:,1));
                Sc(1,1) = mean(S(:,2));
                
                anzahl_kontur = length(edge_index(:,f_coincidential));
                speicher_kontur_oben_3 = zeros(length(edge_index(:,f_coincidential)), 7);
                speicher_kontur_unten_2 = zeros(length(edge_index(:,f_coincidential)), 7);
                coincidential_edge_index = edge_index(:,f_coincidential);
                coincidential_nozzle_pos = obj.position_nozzle(f_coincidential);
                
                perpendicular_edge_index = edge_index(:,f_perpendicular);
                perpendicular_nozzle_pos = obj.position_nozzle(f_perpendicular);
                
                for lauf_kontur = 1:1:anzahl_kontur
                    vect_nozzle_kontur = [(coincidential_edge_index(lauf_kontur) - coincidential_nozzle_pos) (perpendicular_edge_index(lauf_kontur) -perpendicular_nozzle_pos)];
                    vect_nozzle_schwerpunkt = [(Sc(1,f_coincidential) - coincidential_nozzle_pos) (Sc(1,f_perpendicular) -perpendicular_nozzle_pos)];
                    length_vect_nozzle_kontur = sqrt(sum(vect_nozzle_kontur.^2));
                    distance_kontur_nozzle = length_vect_nozzle_kontur;
                    length_vect_nozzle_schwerpunkt = sqrt(sum(vect_nozzle_schwerpunkt.^2));
                    dp_schwerpunkt_nozzle_kontur = dot(vect_nozzle_kontur, vect_nozzle_schwerpunkt);
                    projektionswinkel_kontur = acos(dp_schwerpunkt_nozzle_kontur/(length_vect_nozzle_kontur*length_vect_nozzle_schwerpunkt));
                    distance_projkontur_nozzle = distance_kontur_nozzle * cos(projektionswinkel_kontur);
                    distance_projkontur_kontur = distance_kontur_nozzle * sin(projektionswinkel_kontur);
                    coincidential_proj_kontur = (coincidential_nozzle_pos) + (Sc(1,f_coincidential)-coincidential_nozzle_pos) / length_vect_nozzle_schwerpunkt * distance_projkontur_nozzle;
                    perpendicular_proj_kontur = (perpendicular_nozzle_pos) + (Sc(1,f_perpendicular)-perpendicular_nozzle_pos) / length_vect_nozzle_schwerpunkt * distance_projkontur_nozzle;
                    
                    if (edge_index(lauf_kontur,f_perpendicular)-Sc(1,f_perpendicular))<=0
                        %                 if turnsign * sign(center_nozzle_vec(f_coincidential))*vect_nozzle_schwerpunkt(2) >= 0
                        speicher_kontur_oben_3(lauf_kontur,1) = edge_index(lauf_kontur,1);
                        speicher_kontur_oben_3(lauf_kontur,2) = edge_index(lauf_kontur,2);
                        speicher_kontur_oben_3(lauf_kontur,5) = distance_projkontur_kontur;
                        speicher_kontur_oben_3(lauf_kontur,6) = distance_projkontur_nozzle;
                        speicher_kontur_oben_3(lauf_kontur,7) = projektionswinkel_kontur*180/pi;
                    else
                        speicher_kontur_unten_2(lauf_kontur,1) = edge_index(lauf_kontur,1);
                        speicher_kontur_unten_2(lauf_kontur,2) = edge_index(lauf_kontur,2);
                        speicher_kontur_unten_2(lauf_kontur,5) = distance_projkontur_kontur;
                        speicher_kontur_unten_2(lauf_kontur,6) = distance_projkontur_nozzle;
                        speicher_kontur_unten_2(lauf_kontur,7) = projektionswinkel_kontur*180/pi;
                    end
                end
                speicher_kontur_oben_3_sort = sortrows(speicher_kontur_oben_3,5);
                speicher_kontur_unten_2_sort = sortrows(speicher_kontur_unten_2,5);
                angle_winkelmaxbreite_3 = speicher_kontur_oben_3_sort(anzahl_kontur,7) ;
                angle_winkelmaxbreite_2 = speicher_kontur_unten_2_sort(anzahl_kontur,7) ;
                abstand_nozzle_winkelmaxbreite_3_auf_gerade =  speicher_kontur_oben_3_sort(anzahl_kontur,5) * obj.scale;
                abstand_nozzle_winkelmaxbreite_2_auf_gerade =  speicher_kontur_unten_2_sort(anzahl_kontur,5) * obj.scale;
                %MaKr
                abstand_nozzle_winkelmaxbreite_3_auf_gerade_breite =  speicher_kontur_oben_3_sort(anzahl_kontur,6) * obj.scale;
                abstand_nozzle_winkelmaxbreite_2_auf_gerade_breite =  speicher_kontur_unten_2_sort(anzahl_kontur,6) * obj.scale;
                % MaKr ende
                
                [Kegelwinkel_6_2] = angle_winkelmaxbreite_2;
                [Kegelwinkel_6_3] = angle_winkelmaxbreite_3;
                [Kegelwinkel_6_2_abstand] = abstand_nozzle_winkelmaxbreite_2_auf_gerade;
                [Kegelwinkel_6_3_abstand] = abstand_nozzle_winkelmaxbreite_3_auf_gerade;
                % MaKr
                [Kegelwinkel_6_2_abstand_breite]=abstand_nozzle_winkelmaxbreite_2_auf_gerade_breite;
                [Kegelwinkel_6_3_abstand_breite]=abstand_nozzle_winkelmaxbreite_3_auf_gerade_breite;
                % MaKr ende
                
                plume_angle_max_spray_width1 = [Kegelwinkel_6_2];
                plume_angle_max_spray_width2 = [Kegelwinkel_6_3];
                plume_angle_max_spray_width3 = [Kegelwinkel_6_2_abstand];
                plume_angle_max_spray_width4 = [Kegelwinkel_6_3_abstand];
                % MaKr
                plume_angle_max_spray_width5 = [Kegelwinkel_6_2_abstand_breite];
                plume_angle_max_spray_width6 = [Kegelwinkel_6_3_abstand_breite];
                %MaKr ende
            end
            if isempty(plume_angle_max_spray_width1)
                plume_angle_max_spray_width1 = NaN;
            elseif isempty(plume_angle_max_spray_width2)
                plume_angle_max_spray_width2 = NaN;
            elseif isempty(plume_angle_max_spray_width3)
                plume_angle_max_spray_width3 = NaN;
            elseif isempty(plume_angle_max_spray_width4)
                plume_angle_max_spray_width4 = NaN;
                %MaKr
            elseif isempty(plume_angle_max_spray_width5)
                plume_angle_max_spray_width4 = NaN;
            elseif isempty(plume_angle_max_spray_width5)
                plume_angle_max_spray_width4 = NaN;
                %MaKr ende
            end
            obj.plumeangle_maxspraywidth(imagenumb,seriesnumb) = plume_angle_max_spray_width1+plume_angle_max_spray_width2;
            seriesnumbmax = floor(size(obj.imagenames,2)/obj.Numberofpics);
            if sum(isnan(obj.sprayXmask_crossing(imagenumb,:)))> floor(0.2*seriesnumbmax)
                obj.plumeangle_maxspraywidth(imagenumb,:) = NaN;
            else
                sprayXvector_crossing = mean(obj.sprayXmask_crossing(imagenumb,seriesnumb));
                obj.plumeangle_maxspraywidth(imagenumb,seriesnumb) = obj.plumeangle_maxspraywidth(imagenumb,seriesnumb)*sprayXvector_crossing;
            end
            
            %             storage_plume_angle_max_spray_width1(imagenumb,seriesnumb) = plume_angle_max_spray_width1;
            %             storage_plume_angle_max_spray_width2(imagenumb,seriesnumb) = plume_angle_max_spray_width2;
            %             storage_plume_angle_max_spray_width3(imagenumb,seriesnumb) = plume_angle_max_spray_width3;
            %             storage_plume_angle_max_spray_width4(imagenumb,seriesnumb) = plume_angle_max_spray_width4;
            %             %MaKr
            %             storage_plume_angle_max_spray_width5(imagenumb,seriesnumb) = plume_angle_max_spray_width5;
            %             storage_plume_angle_max_spray_width6(imagenumb,seriesnumb) = plume_angle_max_spray_width6;
        end
        function obj = eval_plumeangle_circle(obj,imagenr,picture_edge)
            circle_radius = 5;
            center_nozzle_vec = [size(obj.bw_mask_spray,1)/2 - obj.position_nozzle(1,1),size(obj.bw_mask_spray,2)/2 - obj.position_nozzle(1,2)];
            seriesnumb = ceil(imagenr/obj.Numberofpics);
            imagenumb = imagenr-(seriesnumb-1)*obj.Numberofpics;
            [row, col] = find(picture_edge); %row = y, col = x
            edge_index = [row, col];
            f_coincidential = find(abs(center_nozzle_vec)==max(abs(center_nozzle_vec)));
            if isempty(row)
                obj.plumeangle_circle(imagenumb,seriesnumb) = NaN;
            else
                distance = round((( obj.position_nozzle(1,1)-row).^2+( obj.position_nozzle(1,2)-col).^2).^0.5);
                d = find((round(circle_radius/obj.scale) - 1) <= distance & distance <= (round(circle_radius/obj.scale) + 1));
                if isempty(d)
                    obj.plumeangle_circle(imagenumb,seriesnumb) = NaN;
                else
                    y1_kreis = max(row(d));
                    y2_kreis = min(row(d));
                    x1_kreis = max(col(d));
                    x2_kreis = min(col(d));
                    %                 x1_kreis = round( -[round((circle_radius/obj.scale)^2) - (obj.position_nozzle(1,1)-y1_kreis)^2].^0.5 + obj.position_nozzle(1,2));
                    %                 x2_kreis = round( -[round((circle_radius/obj.scale)^2) - (obj.position_nozzle(1,1)-y2_kreis)^2].^0.5 + obj.position_nozzle(1,2));
                    vect1_kreis = [ (obj.position_nozzle(1,1)-x1_kreis) (obj.position_nozzle(1,2)-y1_kreis)];
                    vect2_kreis = [ (obj.position_nozzle(1,1)-x2_kreis) (obj.position_nozzle(1,2)-y2_kreis)];
                    obj.plumeangle_circle(imagenumb,seriesnumb) = acos(dot(vect1_kreis, vect2_kreis)/(sqrt(sum(vect1_kreis.^2))*sqrt(sum(vect2_kreis.^2))))*180/pi;
                end
            end
        end
        function obj = eval_spraydirection(obj,imagenr,picture_edge,picture_bwarea,picture_full_adjust_imcom)
            center_nozzle_vec = [size(obj.bw_mask_spray,1)/2 - obj.position_nozzle(1,1),size(obj.bw_mask_spray,2)/2 - obj.position_nozzle(1,2)];
            seriesnumb = ceil(imagenr/obj.Numberofpics);
            imagenumb = imagenr-(seriesnumb-1)*obj.Numberofpics;
            [row, col] = find(picture_edge); %row = y, col = x
            edge_index = [row, col];
            f_coincidential = find(abs(center_nozzle_vec)==max(abs(center_nozzle_vec)));
            if isempty(row)
                obj.spraydirection(imagenumb,seriesnumb) = NaN;
            else
                picture_bw = bwlabel(picture_bwarea);
                S = (struct2vektorWC(regionprops(picture_bw, picture_full_adjust_imcom,'WeightedCentroid')));
                Sc(1,2) = mean(S(:,1));
                Sc(1,1) = mean(S(:,2));
                vect_strahlrichtung_referenz = [ 0 0 ]; % original: vect_strahlrichtung_referenz = [ ( SIT_parameters{9,2} - SIT_parameters{9,2} - 20 ) ( SIT_parameters{9,3} - SIT_parameters{9,3} ) ];
                vect_strahlrichtung_referenz(f_coincidential) = sign(center_nozzle_vec(f_coincidential)) * 20;
                vect_strahlrichtung_winkel = [ ( Sc(1,1) - obj.position_nozzle(1)) ( Sc(1,2) - obj.position_nozzle(2)) ];
                length_strahlrichtung_referenz = sqrt(sum(vect_strahlrichtung_referenz.^2));
                length_strahlrichtung_winkel = sqrt(sum(vect_strahlrichtung_winkel.^2));
                dp_strahlrichtung = dot(vect_strahlrichtung_winkel, vect_strahlrichtung_referenz);
                seriesnumbmax = floor(size(obj.imagenames,2)/obj.Numberofpics);
                obj.spraydirection(imagenumb,seriesnumb) = acos(dp_strahlrichtung/(length_strahlrichtung_referenz * length_strahlrichtung_winkel)) * 180 / pi;
                if sum(isnan(obj.sprayXmask_crossing(imagenumb,:)))> floor(0.2*seriesnumbmax)
                    obj.spraydirection(imagenumb,:) = NaN;
                else
                    sprayXvector_crossing = mean(obj.sprayXmask_crossing(imagenumb,seriesnumb));
                    obj.spraydirection(imagenumb,seriesnumb) = obj.spraydirection(imagenumb,seriesnumb)*sprayXvector_crossing;
                end
            end
        end
        function obj = eval_sprayarea(obj,imagenr,picture_edge,picture_bwarea)
            %             center_nozzle_vec = [size(obj.bw_mask_spray,1)/2 - obj.position_nozzle(1,1),size(obj.bw_mask_spray,2)/2 - obj.position_nozzle(1,2)];
            seriesnumb = ceil(imagenr/obj.Numberofpics);
            imagenumb = imagenr-(seriesnumb-1)*obj.Numberofpics;
            [row, col] = find(picture_edge); %row = y, col = x
            %             edge_index = [row, col];
            %             f_coincidential = find(abs(center_nozzle_vec)==max(abs(center_nozzle_vec)));
            if isempty(row)
                obj.sprayarea(imagenumb,seriesnumb) = NaN;
            else
                obj.sprayarea(imagenumb,seriesnumb) = length(find(~picture_bwarea == 0)) * obj.scale * obj.scale;
                seriesnumbmax = floor(size(obj.imagenames,2)/obj.Numberofpics);
                if sum(isnan(obj.sprayXmask_crossing(imagenumb,:)))> floor(0.2*seriesnumbmax)
                    obj.sprayarea(imagenumb,:) = NaN;
                else
                    sprayXvector_crossing = mean(obj.sprayXmask_crossing(imagenumb,seriesnumb));
                    obj.sprayarea(imagenumb,seriesnumb) = obj.sprayarea(imagenumb,seriesnumb)*sprayXvector_crossing;
                end
                
            end
        end
        
        % Auxiliary functions
        function SheetName = Choose_SheetNames(obj)
            [indx,tf] = listdlg('PromptString','Select Sheetname Values:','SelectionMode','multiple','ListString',obj.Docu_Names);
            for k = 1 : length(indx)
                %                 f0 = strfind(obj.Docu_Names,f{k});
                %                 f0 = find(not(cellfun('isempty', f0)));
                if k == 1
                    SheetName = obj.Docu_Values{indx};
                else
                    SheetName = [SheetName,'_',obj.Docu_Values{indx}];
                end
                
            end
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
        
        % Export functions
        function picture_overlay = Export_ImageEdge(obj,imagenr,picture_edge,picture_normed,switchexport)
            bw_mask_spray_edge = edge(obj.bw_mask_spray, 'canny');
            picture_mask_edge =   picture_edge;
            picture_mask_nozzle = zeros(size(picture_mask_edge));
            picture_mask_nozzle(round(obj.position_nozzle(1,1)),round(obj.position_nozzle(1,2))) = 1;
            picture_mask_B = picture_normed.*(1-picture_mask_edge).*(1-bw_mask_spray_edge).*(1-picture_mask_nozzle)+(1-picture_normed).*bw_mask_spray_edge;
            picture_mask_G = picture_normed.*(1-picture_mask_edge).*(1-bw_mask_spray_edge).*(1-picture_mask_nozzle)+(1-picture_normed).*picture_mask_nozzle;
            picture_mask_R = picture_normed.*(1-picture_mask_edge).*(1-bw_mask_spray_edge).*(1-picture_mask_nozzle)+(1-picture_normed).*picture_mask_edge;
            picture_overlay = uint8(cat(3,picture_mask_R,picture_mask_G,picture_mask_B)*255);
            
            %                 folder_export = strrep(obj.folder,'18-01 PN Glasdüse Soligor','18-01 PN Glasdüse Soligor Export Edge');
            if switchexport == 1
                folder_export = strrep(obj.folder,'18-06 PN Glass Nozzle Min Max jwp','18-06 PN Glass Nozzle Min Max jwp Export Edge');
                folder_export = [obj.current_path{1},'\',folder_export{1}];
                %             folder_export = ['F:','\',folder_export{1}];
                if ~isfolder(folder_export)
                    mkdir(folder_export);
                end
                if ~isfile([folder_export,'\',obj.imagenames{imagenr}])
                    imwrite(picture_overlay,[folder_export,'\',obj.imagenames{imagenr}]);
                else
                    disp(['Warning: ',folder_export,'\',obj.imagenames{imagenr},' does already exist. Image was not written']);
                end
            end
        end
        function obj = Export_Excel(varargin)
            obj = varargin{1};
            if nargin > 1
                master_path = varargin{2};
            else
                [FileName,PathName,FilterIndex] = uigetfile;
                master_path = [PathName,FileName];
            end
            obj = check_path(obj);
            destination = [pwd,'\Excel_Export'];
            SheetName = Choose_SheetNames(obj);
            if ~isfolder(destination)
                mkdir(destination);
            end
            property_names = properties(obj);
            for i = 1 : length(property_names)
                if ~isfile([destination,'\',property_names{i},'.xlsx'])
                    copyfile(master_path,[destination,'\',property_names{i},'.xlsx']);
                end
                excel = actxserver('excel.application');
                wkbk = excel.Workbooks.Open([destination,'\',property_names{i},'.xlsx']);
                wksheet = wkbk.Worksheets.Item('Master');
                wksheet.Copy(wksheet);
                newSheet=wkbk.Worksheets.Item('Master (2)');
                newSheet.Name=SheetName;
                wkbk.Save
                excel.Quit
                excel.delete
                Data_Cell = transpose(obj.(char(property_names{i})));
                xlswrite([destination,'\',property_names{i},'.xlsx'],obj.time,SheetName,'B18');
                xlswrite([destination,'\',property_names{i},'.xlsx'],Data_Cell,SheetName,'B20');
                
            end
            
        end
        function Data = Export_Data(obj)
            Data = struct(obj);
            %             propertylist = properties(obj);
            %             for i = 1 : length(propertylist)
            %                 Data.(char(propertylist{i})) = obj.(char(propertylist{i}));
            %             end
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
        function time = get.time(obj)
            time = 0:10^6 / (1000*obj.frequency):(10^6 / (1000*obj.frequency)*obj.Numberofpics-10^6 / (1000*obj.frequency));
        end
        
        % Figure functions
        function argout = timeplot(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            addRequired(p,'obj');
            addRequired(p,'cfh');
            addRequired(p,'clh');
            addParameter(p,'Color','b');
            addParameter(p,'Legend',{'default'});
            addParameter(p,'YData','default');
            addParameter(p,'LineStyle','-');
            addParameter(p,'XTime',0);
            addParameter(p,'mode','default');
            addParameter(p,'MarkerStyle','x');
            
            
            parse(p,varargin{:});
            obj = p.Results.obj;
            cfh = p.Results.cfh;
            clh = p.Results.clh;
            Color = p.Results.Color;
            Legend = p.Results.Legend;
            YData = p.Results.YData;
            XTime = p.Results.XTime;
            LineStyle = p.Results.LineStyle;
            mode = p.Results.mode;
            MarkerStyle = p.Results.MarkerStyle;
            
            Legendstring = '';
            
            if strcmp(mode,'default')==1
                output_ok = 0;
                while output_ok == 0
                mode = questdlg('Please choose how you want to plot your data','Input', 'mean','median', 'single','mean');
                if ~isempty(mode)
                    output_ok = 1;
                end
                end
            end
            
            if strcmp(Legend{1},'default')==1
                [indx,tf] = listdlg('PromptString','Select Legendentries:','SelectionMode','multiple','ListString',obj.Docu_Names);
                Legend = obj.Docu_Names(indx);
            end
            
            if strcmp(YData,'default')==1
                listofproperties = properties(obj);
                [indx,tf] = listdlg('PromptString','Select YData:','SelectionMode','single','ListString',listofproperties);
                YData = listofproperties{indx};
            end
            
            if sum(XTime)==0
                listoftime = sprintfc('%d',obj.time);
                [indx,tf] = listdlg('PromptString','Select relevant timesteps for XData','SelectionMode','multiple','ListString',listoftime);
                XTime = indx;
            end
            
            for i = 1:length(Legend)
                for j = 1 : length(obj.Docu_Names)
                    if strcmp(Legend{i},obj.Docu_Names{j})
                        Legendstring = [Legendstring,' ',obj.Docu_Values{j}];
                    end
                end
            end
            
            set(0, 'currentfigure', cfh);  %# for figures
            grid('on');
            hold on
            % set(f, 'currentaxes', axs);  %# for axes with handle axs on figure f
            default_cmap = parula(size(obj.(char(YData)),2));
            if strcmp(LineStyle,'NoLine')
%                 scatter(nanmean(X), nanmean(Y), 'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',Color,'LineWidth',1);
                if strcmp(mode,'mean')
                    scatter(obj.time(1,XTime)', mean(obj.(char(YData))(XTime,:),2), MarkerStyle,'MarkerEdgeColor', Color,'LineWidth', 1);
                elseif strcmp(mode,'median')
                    scatter(obj.time(1,XTime)', median(obj.(char(YData))(XTime,:),2,'omitnan'),MarkerStyle, 'MarkerEdgeColor', Color,'LineWidth', 1);
                elseif strcmp(mode,'single')
                    for i = 1 : size(obj.(char(YData)),2)
                        scatter(obj.time(1,XTime)', obj.(char(YData))(XTime,i), MarkerStyle,'MarkerEdgeColor', default_cmap(i,:),'LineWidth', 1);
                    end
                end
            else
                if strcmp(mode,'mean')
                    plot(obj.time(1,XTime)', mean(obj.(char(YData))(XTime,:),2), 'Color', Color, 'LineStyle', LineStyle, 'LineWidth', 1);
                elseif strcmp(mode,'median')
                    plot(obj.time(1,XTime)', median(obj.(char(YData))(XTime,:),2,'omitnan'), 'Color', Color, 'LineStyle', LineStyle, 'LineWidth', 1);
                elseif strcmp(mode,'single')
                    for i = 1 : size(obj.(char(YData)),2)
                        plot(obj.time(1,XTime)', obj.(char(YData))(XTime,i), 'Color', default_cmap(i,:), 'LineStyle', LineStyle, 'LineWidth', 1);
                    end
                end
            end
            set(0, 'currentfigure', clh);  %# for figures
            hold on
            hLegend = findobj(gcf, 'Type', 'Legend');
            % scatter
            if strcmp(LineStyle,'NoLine')
                if isempty(hLegend)
                    if strcmp(mode,'mean') || strcmp(mode,'median')
                        scatter(NaN,NaN, MarkerStyle,'MarkerEdgeColor', Color,'LineWidth', 1);
                        legend(Legendstring, 'Location', 'Best', 'Orientation', 'Vertical');
                    elseif strcmp(mode,'single')
                        Legendcell = {};
                        for i = 1 : size(obj.(char(YData)),2)
                            scatter(0,0, MarkerStyle,'MarkerEdgeColor', default_cmap(i,:),'LineWidth', 1);
                            Legendcell = [Legendcell,[Legendstring,' ShotNr ',num2str(i)]];
                        end
                        legend(Legendcell, 'Location', 'Best', 'Orientation', 'Vertical');
                    end
                    
                else
                    past_Legend = hLegend.String;
                    if strcmp(mode,'mean') || strcmp(mode,'median')
                        scatter(NaN,NaN, MarkerStyle,'MarkerEdgeColor', Color,'LineWidth', 1);
                        legend([past_Legend,Legendstring], 'Location', 'Best', 'Orientation', 'Vertical');
                    elseif strcmp(mode,'single')
                        Legendcell = past_Legend;
                        for i = 1 : size(obj.(char(YData)),2)
                            scatter(0,0, MarkerStyle,'MarkerEdgeColor', default_cmap(i,:),'LineWidth', 1);
                            Legendcell = [Legendcell,[Legendstring,' ShotNr ',num2str(i)]];
                        end
                        legend(Legendcell, 'Location', 'Best', 'Orientation', 'Vertical');
                    end
                end
            % plot
            else
                if isempty(hLegend)
                    if strcmp(mode,'mean') || strcmp(mode,'median')
                        plot(0, 0,  'Color', Color, 'LineStyle', LineStyle, 'LineWidth', 1);
                        legend(Legendstring, 'Location', 'Best', 'Orientation', 'Vertical');
                    elseif strcmp(mode,'single')
                        Legendcell = {};
                        for i = 1 : size(obj.(char(YData)),2)
                            plot(0, 0, 'Color', default_cmap(i,:), 'LineStyle', LineStyle, 'LineWidth', 1);
                            Legendcell = [Legendcell,[Legendstring,' ShotNr ',num2str(i)]];
                        end
                        legend(Legendcell, 'Location', 'Best', 'Orientation', 'Vertical');
                    end
                    
                else
                    past_Legend = hLegend.String;
                    if strcmp(mode,'mean') || strcmp(mode,'median')
                        plot(0, 0,  'Color', Color, 'LineStyle', LineStyle, 'LineWidth', 1);
                        legend([past_Legend,Legendstring], 'Location', 'Best', 'Orientation', 'Vertical');
                    elseif strcmp(mode,'single')
                        Legendcell = past_Legend;
                        for i = 1 : size(obj.(char(YData)),2)
                            plot(0, 0, 'Color', default_cmap(i,:), 'LineStyle', LineStyle, 'LineWidth', 1);
                            Legendcell = [Legendcell,[Legendstring,' ShotNr ',num2str(i)]];
                        end
                        legend(Legendcell, 'Location', 'Best', 'Orientation', 'Vertical');
                    end
                end
            end
            argout{1} = 'Legend';
            argout{2} = Legend;
            argout{3} = 'YData';
            argout{4} = YData;
            argout{5} = 'XTime';
            argout{6} = XTime;
            argout{7} = 'mode';
            argout{8} = mode;
        end
        function argout = scatterplot(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            addRequired(p,'obj');
            addRequired(p,'cfh');
            addRequired(p,'clh');
            addParameter(p,'Color','b');
            addParameter(p,'Legend',{'default'});
            addParameter(p,'YData','default');
            addParameter(p,'XData','default');
            addParameter(p,'YTime',0);
            addParameter(p,'XTime',0);
            addParameter(p,'ErrBar','off');
            addParameter(p,'mkr','o');
            addParameter(p,'NrofOP',1);
            addParameter(p,'CurrentNrofOP',1);
            
            parse(p,varargin{:});
            obj = p.Results.obj;
            cfh = p.Results.cfh;
            clh = p.Results.clh;
            Color = p.Results.Color;
            Legend = p.Results.Legend;
            YData = p.Results.YData;
            XData = p.Results.XData;
            YTime = p.Results.YTime;
            XTime = p.Results.XTime;
            ErrBar = p.Results.ErrBar;
            mkr = p.Results.mkr;
            NrofOP = p.Results.NrofOP;
            CurrentNrofOP = p.Results.CurrentNrofOP;
            
            Legendstring = '';
            if strcmp(Legend{1},'default')==1
                [indx,tf] = listdlg('PromptString','Select Legendentries:','SelectionMode','multiple','ListString',obj.Docu_Names);
                Legend = obj.Docu_Names(indx);
            end
            
            if strcmp(YData,'default')==1
                listofproperties = properties(obj);
                [indx,tf] = listdlg('PromptString','Select YData','SelectionMode','single','ListString',listofproperties);
                YData = listofproperties{indx};
            end
            
            if sum(YTime)==0
                listoftime = sprintfc('%d',obj.time);
                [indx,tf] = listdlg('PromptString','Select relevant timesteps for YData','SelectionMode','multiple','ListString',listoftime);
                YTime = indx;
            end
            
            if strcmp(XData,'default')==1
                listofproperties = properties(obj);
                [indx,tf] = listdlg('PromptString','Select XData','SelectionMode','single','ListString',listofproperties);
                XData = listofproperties{indx};
            end
            if sum(XTime)==0
                listoftime = sprintfc('%d',obj.time);
                [indx,tf] = listdlg('PromptString','Select relevant timesteps for XData','SelectionMode','multiple','ListString',listoftime);
                XTime = indx;
            end
            
            for i = 1:length(Legend)
                for j = 1 : length(obj.Docu_Names)
                    if strcmp(Legend{i},obj.Docu_Names{j})
                        Legendstring = [Legendstring,' ',obj.Docu_Values{j}];
                    end
                end
            end
            
            set(0, 'currentfigure', cfh);  %# for figures
            grid('on');
            hold on
            % set(f, 'currentaxes', axs);  %# for axes with handle axs on figure f
            Y = obj.(char(YData))(YTime,:);
            Y = reshape(Y,[size(Y,1)*size(Y,2),1]);
            X = obj.(char(XData))(XTime,:);
            X = reshape(X,[size(X,1)*size(X,2),1]);
            scatter(nanmean(X), nanmean(Y), mkr,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',Color,'LineWidth',1);
            if strcmp(ErrBar,'on')
                yneg = nanstd(Y);
                ypos = nanstd(Y);
                xneg = nanstd(X);
                xpos = nanstd(X);
                errorbar(nanmean(X), nanmean(Y),yneg,ypos,xneg,xpos,'Color',[0 0 0])
            end
            set(0, 'currentfigure', clh);  %# for figures
            hold on
            hLegend = findobj(gcf, 'Type', 'Legend');
            if isempty(hLegend)
                scatter(0, 0, mkr,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',Color,'LineWidth',1);
                legend(Legendstring, 'Location', 'Best', 'Orientation', 'Vertical');
                legend('boxoff')
            else
                past_Legend = hLegend.String;
                scatter(0, 0, mkr,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',Color,'LineWidth',1);
                legend([past_Legend,Legendstring], 'Location', 'Best', 'Orientation', 'Vertical');
                legend('boxoff')
            end
            
            set(0, 'currentfigure', cfh);  %# for figures
            if NrofOP == CurrentNrofOP
                test =1;
            end
            
            argout{1} = 'Legend';
            argout{2} = Legend;
            argout{3} = 'YData';
            argout{4} = YData;
            argout{5} = 'XData';
            argout{6} = XData;
            argout{7} = 'YTime';
            argout{8} = YTime;
            argout{9} = 'XTime';
            argout{10} = XTime;
            argout{11} = 'ErrBar';
            argout{12} = ErrBar;
            argout{13} = 'mkr';
            argout{14} = mkr;
        end
        function argout = box_plot(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            addRequired(p,'obj');
            addRequired(p,'cfh');
            addRequired(p,'clh');
            addParameter(p,'Color','b');
            addParameter(p,'BoxLegend',{'default'});
            addParameter(p,'YData','default');
            addParameter(p,'XData','default');
            addParameter(p,'YTime',0);
            addParameter(p,'XTime',0);
            addParameter(p,'NrofOP',1);
            addParameter(p,'CurrentNrofOP',1);
            addParameter(p,'boxdata',[]);
            addParameter(p,'boxlabel',[]);
            
            
            parse(p,varargin{:});
            obj = p.Results.obj;
            cfh = p.Results.cfh;
            clh = p.Results.clh;
            Color = p.Results.Color;
            BoxLegend = p.Results.BoxLegend;
            YData = p.Results.YData;
            YTime = p.Results.YTime;
            NrofOP = p.Results.NrofOP;
            CurrentNrofOP = p.Results.CurrentNrofOP;
            boxdata = p.Results.boxdata;
            boxlabel = p.Results.boxlabel;
            
            if strcmp(BoxLegend{1},'default')==1
                [indx,tf] = listdlg('PromptString','Select Legendentries:','SelectionMode','multiple','ListString',obj.Docu_Names);
                BoxLegend = obj.Docu_Names(indx);
            end
            
            if strcmp(YData,'default')==1
                listofproperties = properties(obj);
                [indx,tf] = listdlg('PromptString','Select YData','SelectionMode','single','ListString',listofproperties);
                YData = listofproperties{indx};
            end
            
            if sum(YTime)==0
                listoftime = sprintfc('%d',obj.time);
                [indx,tf] = listdlg('PromptString','Select relevant timesteps for YData','SelectionMode','multiple','ListString',listoftime);
                YTime = indx;
            end
            
            groups = cell(size(BoxLegend));
            for i = 1:length(BoxLegend)
                for j = 1 : length(obj.Docu_Names)
                    if strcmp(BoxLegend{i},obj.Docu_Names{j})                
                        if i == 1
                            groups(i) = obj.Docu_Values(j);
                        else
                            groups(i) = obj.Docu_Values(j);
                        end
                    end
                end
            end
            
            %             current_axes = findall(cfh, 'type', 'axes');
            Y = obj.(char(YData))(YTime,:);
            Y = reshape(Y,[size(Y,1)*size(Y,2),1]);
            Y_filtered = Y(~isnan(Y));
            
            %             if isempty(current_axes.UserData)
            if CurrentNrofOP == 1 && isempty(boxdata)
                boxdata = {Y_filtered};
                boxlabel = {groups};
            else
                boxdata = [boxdata,{Y_filtered}];
                boxlabel = [boxlabel,{groups}];
            end
            
            
            if CurrentNrofOP == NrofOP
                set(0, 'currentfigure', cfh);  %# for figures
                grid('on');
                data = [];
                group = {};
                for i = 1 : size(boxdata,2)
                    if i == 1
                        data = boxdata{1,i}';
                        group = repmat(boxlabel{1,i}',1,numel(boxdata{1,i}'));
                    else
                        data = [data,boxdata{1,i}'];
                        group = [group,repmat(boxlabel{1,i}',1,numel(boxdata{1,i}'))];
                    end
                end
                grouping = cell(size(group,1),1);
                for i = 1:size(group,1)
                    grouping{i,1} = group(i,:)';
                end
                prompt = BoxLegend;
                title = 'Rearrange Labels';
                dims = ones(size(BoxLegend));
                definput = 1:1:length(BoxLegend);
                definput = cellstr(num2str(definput'));
%                 definput = {'20','hsv'};
                answer = inputdlg(prompt,title,dims,definput);
                sorting = cellfun(@str2num,answer);
                grouping_sorted = grouping;
                for k = 1 : length(sorting)
                    grouping_sorted(sorting(k)) = grouping(k);
                end
                set(0, 'currentfigure', cfh);
                boxplot(data',grouping_sorted,'LabelVerbosity','minor');
            end
            
            argout{1} = 'YData';
            argout{2} = YData;
            argout{3} = 'YTime';
            argout{4} = YTime;
            argout{5} = 'boxdata';
            argout{6} = boxdata;
            argout{7} = 'boxlabel';
            argout{8} = boxlabel;
            argout{9} = 'BoxLegend';
            argout{10} = BoxLegend;

            
        end
        function argout = histogram_plot(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            addRequired(p,'obj');
            addRequired(p,'cfh');
            addRequired(p,'clh');
            addParameter(p,'Color','b');
            addParameter(p,'Legend',{'default'});
            addParameter(p,'Data','default');
            addParameter(p,'Time',0);
            addParameter(p,'XTime',0);
            addParameter(p,'NrofOP',1);
            addParameter(p,'CurrentNrofOP',1);
            addParameter(p,'NrofBins',[]);
            addParameter(p,'fittype','normal');
            addParameter(p,'multiTile','yes');
            addParameter(p,'BinWidth',2);
            

            parse(p,varargin{:});
            obj = p.Results.obj;
            cfh = p.Results.cfh;
            clh = p.Results.clh;
            Color = p.Results.Color;
            Data = p.Results.Data;
            Time = p.Results.Time;
            NrofOP = p.Results.NrofOP;
            CurrentNrofOP = p.Results.CurrentNrofOP;
            Legend = p.Results.Legend;
            NrofBins = p.Results.NrofBins;
            fittype= p.Results.fittype;
            multiTile = p.Results.multiTile;
            BinWidth = p.Results.BinWidth;
            
            X_Tiles = 5;
             Legendstring = '';
            if strcmp(Legend{1},'default')==1
                [indx,tf] = listdlg('PromptString','Select Legendentries:','SelectionMode','multiple','ListString',obj.Docu_Names);
                Legend = obj.Docu_Names(indx);
            end
            
            if strcmp(Data,'default')==1
                listofproperties = properties(obj);
                [indx,tf] = listdlg('PromptString','Select YData','SelectionMode','single','ListString',listofproperties);
                Data = listofproperties{indx};
            end
            
            if sum(Time)==0
                listoftime = sprintfc('%d',obj.time);
                [indx,tf] = listdlg('PromptString','Select relevant timesteps for YData','SelectionMode','multiple','ListString',listoftime);
                Time = indx;
            end
            
            for i = 1:length(Legend)
                for j = 1 : length(obj.Docu_Names)
                    if strcmp(Legend{i},obj.Docu_Names{j})
                        Legendstring = [Legendstring,' ',obj.Docu_Values{j}];
                    end
                end
            end
            set(0, 'currentfigure', cfh);  %# for figures
            hold on
            Y = obj.(char(Data))(Time,:);
            Y = reshape(Y,[size(Y,1)*size(Y,2),1]);
            if strcmp(multiTile,'yes')
                ax = subplot(ceil(NrofOP/X_Tiles),X_Tiles,CurrentNrofOP);
                set(cfh, 'currentaxes', ax);  %# for axes with handle axs on figure f
            end
            Y = Y(~isnan(Y));
            if isempty(NrofBins)
                NrofBins = ceil((max(Y)-min(Y))/BinWidth);
            end
            
            [muhat,sigma] = normfit(Y);
            % construct histogram
            [F,X] = hist(Y,NrofBins);
            % get ready to plot as probability histogram
            F = F/trapz(X,F);
            b = bar(X,F);
            b.FaceColor = Color;
            grid('on');
            hold on;
            % use muhat and sigma to construct pdf
            x = muhat-3*sigma:0.01:muhat+3*sigma;
            % plot PDF over histogram
            y = normpdf(x,muhat,sigma);
            plot(x,y,'r','linewidth',1);
%            h = histfit(Y,NrofBins,fittype);
%            h(1).FaceColor = Color;
           leg = legend(b,Legendstring,'Location','NorthEast','Orientation','horizontal');
           legend('boxoff')
           rectfg =  get(ax,'Position');
           rect = get(leg,'Position');
           rect = [rectfg(1)+rectfg(3)-rect(3),rectfg(2)+rectfg(4)-rect(4),rect(3),rect(4)];
           set(leg,'Position',rect);
           argout{1} = 'Legend';
           argout{2} = Legend;
           argout{3} = 'Data';
           argout{4} = Data;
           argout{5} = 'Time';
           argout{6} = Time;
        end
        function full_video(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            addRequired(p,'obj');
            addParameter(p,'size',[15,23]);
            addParameter(p,'unittype','centimeters');
            addParameter(p,'seriesnumbstart',1);
            addParameter(p,'seriesnumbend',3);
            addParameter(p,'Name',[strrep(strrep(datestr(datetime('now')),' ','_'),':',''),'_newVideo']);
            addParameter(p,'FontSize',10);
            addParameter(p,'FontName','Calibri');
            parse(p,varargin{:});
            obj = p.Results.obj;
            size = p.Results.size;
            unittype = p.Results.unittype;
            seriesnumbstart = p.Results.seriesnumbstart;
            seriesnumbend = p.Results.seriesnumbend;
            Name = p.Results.Name;
            FontSize = p.Results.FontSize;
            FontName = p.Results.FontName;
            Colors = {'b';'g';'r'};
            plotheight = 0.28;
            plotheightactual = 0.2;
            plotlengthactual = 0.4;
            plotdif = 0.02;
            einzug = 0;
            einzug2 = 0.57;
            v = VideoWriter([obj.current_path{1,1},'\',Name,'.avi']);
            v.FrameRate = 10;
            open(v);
            %             einzug = (1-3*plotheight-2*plotdif)/2;
            back_mean_normed = eval_background(obj,floor(length(obj.time)/2));
            picture_normed_crop =eval_imageadjust(obj,floor(length(obj.time)/2),back_mean_normed);
            [cropping,rect] = imcrop(picture_normed_crop);
                               
            cfh = figure('Name','Current Figure');
            set(cfh,'Units', unittype, 'OuterPosition', [0 0 size(2) size(1)],'color','w');
            ax1 = axes('Position', [0.05 0.05 0.45 0.95]);
            ax2 = axes('Position', [einzug2 1-einzug-plotheight plotlengthactual plotheightactual]);
            ax3 = axes('Position', [einzug2 1-einzug-2*plotheight-plotdif plotlengthactual plotheightactual]);
            ax4 = axes('Position', [einzug2 1-einzug-3*plotheight-2*plotdif plotlengthactual plotheightactual]);
            
            for j = seriesnumbstart : seriesnumbend
                startpoint = 1;
                for i = 1 : length(obj.time)
                    if startpoint == 1
                        if obj.sprayXmask_crossing(i,j) == 1
%                             cfh = figure('Name','Current Figure');
%                             set(cfh,'Units', unittype, 'OuterPosition', [0 0 size(2) size(1)]);
                            imagenr = i+(j-1)*obj.Numberofpics;
                            back_mean_normed = eval_background(obj,imagenr);
                            [picture_normed,picture_full_adjust_imcom] = eval_imageadjust(obj,imagenr,back_mean_normed);
                            [picture_edge,picture_bwarea] = eval_binarize(obj,picture_full_adjust_imcom);
                            picture_normed(picture_bwarea==0)=1;
%                                             picture_overlay = Export_ImageEdge(obj,imagenr,picture_edge,picture_normed,0);
                            picture_normed_cut = imcrop(picture_normed,rect);
%                             picture_normed_cut_rot = imrotate(picture_normed_cut,45);
                            %% Image
%                             ax1 = axes('Position', [0.05 0.05 0.45 0.95]);
                            set(cfh, 'currentaxes', ax1);  %# for axes with handle axs on figure f
                          imshow(picture_normed_cut);
                            %% Plot 1
%                             ax2 = axes('Position', [einzug2 1-einzug-plotheight plotlengthactual plotheightactual]);
                            set(cfh, 'currentaxes', ax2);
                            plot(obj.time(1:i)', obj.total_penetration(1:i,j), 'Color',Colors{j}, 'LineStyle', '-', 'LineWidth', 1);
                            
                            xlabel('Time in ms', 'FontSize', FontSize, 'FontName', FontName); %Schriftgröße und Art Achsenttitel x
                            ylabel('Total Penetration in mm', 'FontSize', FontSize, 'FontName', FontName); %Schriftgröße und Art Achsenttitel y
                            
                            %  set('XTick', Xlim); % Setzt die Anzahl der Teilstriche der x-Achse
                            mult = ~isnan(obj.total_penetration(:,seriesnumbstart : seriesnumbend));
                            set(ax2,'XLim', [0,max(max(obj.time'.*mult))]); % Setzt die Anzahl der Teilstriche der x-Achse
                            hold on
                            % set('YTick', Ylim); % Setzt die Anzahl der Teilstriche der y-Achse
                            set(ax2,'YLim', [0,ceil(max(max(obj.total_penetration(:,seriesnumbstart : seriesnumbend)))/10)*10]); % Setzt die Anzahl der Teilstriche der y-Achse
                            %% Plot 2
%                             ax3 = axes('Position', [einzug2 1-einzug-2*plotheight-plotdif plotlengthactual plotheightactual]);
                            set(cfh, 'currentaxes', ax3);
                            plot(obj.time(1:i)', obj.plumeangle_maxspraywidth(1:i,j), 'Color', Colors{j}, 'LineStyle', '-', 'LineWidth', 1);
                            xlabel('Time in ms', 'FontSize', FontSize, 'FontName', FontName); %Schriftgröße und Art Achsenttitel x
                            ylabel('Split Angle in °', 'FontSize', FontSize, 'FontName', FontName); %Schriftgröße und Art Achsenttitel y
                            
                            %  set('XTick', Xlim); % Setzt die Anzahl der Teilstriche der x-Achse
                            mult = ~isnan(obj.total_penetration(:,seriesnumbstart : seriesnumbend));
                            set(ax3,'XLim', [0,max(max(obj.time'.*mult))]); % Setzt die Anzahl der Teilstriche der x-Achse
                            hold on
                            % set('YTick', Ylim); % Setzt die Anzahl der Teilstriche der y-Achse
                            set(ax3,'YLim', [0,ceil(max(max(obj.plumeangle_maxspraywidth(:,seriesnumbstart : seriesnumbend)))/10)*10]); % Setzt die Anzahl der Teilstriche der y-Achse
                            %% Plot3
%                             ax4 = axes('Position', [einzug2 1-einzug-3*plotheight-2*plotdif plotlengthactual plotheightactual]);
                            set(cfh, 'currentaxes', ax4);
                            plot(obj.time(1:i)', obj.sprayarea(1:i,j)/100, 'Color', Colors{j}, 'LineStyle', '-', 'LineWidth', 1);
                            xlabel('Time in ms', 'FontSize', FontSize, 'FontName', FontName); %Schriftgröße und Art Achsenttitel x
                            ylabel('Spray Area in cm^2', 'FontSize', FontSize, 'FontName', FontName); %Schriftgröße und Art Achsenttitel y
                            
                            %  set('XTick', Xlim); % Setzt die Anzahl der Teilstriche der x-Achse
                            mult = ~isnan(obj.total_penetration(:,seriesnumbstart : seriesnumbend));
                            set(ax4,'XLim', [0,max(max(obj.time'.*mult))]); % Setzt die Anzahl der Teilstriche der x-Achse
                            
                            % set('YTick', Ylim); % Setzt die Anzahl der Teilstriche der y-Achse
                            set(ax4,'YLim', [0,ceil(max(max(obj.sprayarea(:,seriesnumbstart : seriesnumbend)/100))/1)*1]); % Setzt die Anzahl der Teilstriche der y-Achse
                            hold on
                            F = getframe(cfh);
                            writeVideo(v,F);
%                             close(cfh);
                        else
                            startpoint = 0;
                        end
                    end
                end
            end
            close(v);
            close(cfh)
        end
        
        
    end
end


