function [strucoutput] = Parameter_func(strucinput,varargin)
%UNTITLED Summary of this function goes here
% profile: {'HS','DoubleShutter','SingleShutter','Spray','Nozzle'}
%   Detailed explanation goes here
p = inputParser;
p.KeepUnmatched = true;
addRequired(p,'strucinput',@isstruct);
addParameter(p,'oneforall',0);
default_template = struct;
addParameter(p,'template',default_template);
addParameter(p,'Mode','');
addParameter(p,'Current_path','');
parse(p,strucinput,varargin{:});
strucinput = p.Results.strucinput;
oneforall = p.Results.oneforall;
current_template = p.Results.template;
Mode = p.Results.Mode;
Current_path = p.Results.Current_path;
template_exists = 0;
if ~isempty(fieldnames(current_template))
  current_template.Copy_Template = 1;
  template_exists = 1;
end
EvaluationMode_Path = [pwd,'\Evaluation Mode'];
addpath(EvaluationMode_Path);
EvaluationMode_Path = dir(EvaluationMode_Path);
if strcmp(Mode,'')
[indx,tf] = listdlg('PromptString','Select an Evaluation Mode:','SelectionMode','single','ListString',{EvaluationMode_Path.name});
fh = str2func(strrep(EvaluationMode_Path(indx).name,'.m',''));
else
    fh = str2func(strrep(Mode,'.m',''));
end
if strcmp(Mode,'')
Current_Data_Drive = uigetdir('C:\','Select current Data Path');
else
    Current_Data_Drive = Current_path;
end
OPnames = fieldnames(strucinput);



for i=1:length(OPnames)
    
    % Preview of OP
    if  template_exists == 1 && oneforall==1
    else
    winopen([Current_Data_Drive,'\',strucinput.(char(OPnames(i))).folder{1,1}]);
    end
    image =  imread([Current_Data_Drive,'\',strucinput.(char(OPnames(i))).folder{1,1},'\',strucinput.(char(OPnames(i))).imagenames{1,2}]);
    image = double(image);
    image_normed = image./max(max(image));
    image_double = image/2^4;
    image = uint8(image_double);
    
    figure('NumberTitle','off','Name',strucinput.(char(OPnames(i))).imagenames{1,2});
    subplot(1,2,1), imshow(image_normed)
    subplot(1,2,2), imshow(image)
    if template_exists == 1
        if oneforall==1
            loadpara = 'Copy last OP';
        else
            loadpara = questdlg(['Please choose how you want to load your parameters for' strucinput.(char(OPnames(i))).imagenames{1,2} '.'],'Input', 'Manual','Copy last OP','Remove OP','Manual');
        end
    else
        loadpara = questdlg(['Please choose how you want to load your parameters for' strucinput.(char(OPnames(i))).imagenames{1,2} '.'],'Input', 'Manual','Remove OP','Manual');
    end
    close(gcf);
    
    if strcmp(loadpara, 'Remove OP') == 1
                strucinput = rmfield(strucinput,OPnames(i));
    else
   
    if template_exists == 0
        [strucoutput.(char(OPnames(i)))] = fh(strucinput.(char(OPnames(i))));
    else
        if strcmp(loadpara, 'Copy last OP') == 1
            current_template.Copy_Template = 1;
        elseif strcmp(loadpara, 'Manual') == 1
            current_template.Copy_Template = 0;
        end
        [strucoutput.(char(OPnames(i)))] = fh(strucinput.(char(OPnames(i))),'copytemplate',current_template);
    end
    current_template = Export_Data(strucoutput.(char(OPnames(i))));
    template_exists = 1;
    end
end



% OPnames = fieldnames(strucinput);
% GET_Parameters = struct;
% for i=1:length(OPnames)
%     folder = strucinput.(char(OPnames(i))).folder;
%     imagenames = strucinput.(char(OPnames(i))).imagenames;
%     current_path = strucinput.(char(OPnames(i))).current_path;
%     
%     for j = 1 : size(folder,1)
%     current_folder = [current_path{j,1},folder{j,1}];
%     current_imagenames = imagenames(j,:);
%     image = imread([current_folder,'\',current_imagenames{1,1}]);
%     image_orig = double(image);
%     image_normed = pic./max(max(image_orig));
%     image_8bit = uint8(image_normed * 2^8);
%     
%     figure('NumberTitle','off','Name',current_imagenames{1,1});
%     subplot(1,2,1), imshow(image_normed)
%     subplot(1,2,2), imshow(image_8bit)
%     loadpara = questdlg(['Please choose how you want to load your parameters for' imagenames{1} '. Cancel will Skip the OP'],'Input', 'From Workspace', 'From File','Manual','From Workspace');
%     close(gcf);
%     end
% end

