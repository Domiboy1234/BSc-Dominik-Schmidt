function [strucoutput] = Filter_func(strucinput)
%UNTITLED Summary of this function goes here
% profile: {'HS','DoubleShutter','SingleShutter','Spray','Nozzle'}
%   Detailed explanation goes here
EvaluationMode_Path = [pwd,'\Evaluation Mode'];
addpath(EvaluationMode_Path);
EvaluationMode_Path = dir(EvaluationMode_Path);
% [indx,tf] = listdlg('PromptString','Select an Evaluation Mode:','SelectionMode','single','ListString',{EvaluationMode_Path.name});
% fh = str2func(strrep(EvaluationMode_Path(indx).name,'.m',''));
% Current_Data_Drive = uigetdir('C:\','Select current Data Path');
OPnames = fieldnames(strucinput);
% template_exists = 0;

for i=1:length(OPnames)
    if i == 1
        choosemethod = methods(strucinput.(char(OPnames(i))));
        [indx,tf] = listdlg('PromptString','Select a Method for Execution:','SelectionMode','single','ListString',choosemethod);
        fh = str2func(choosemethod{indx});
        [Switch,argout] = fh(strucinput.(char(OPnames(i))));
    else
        [Switch,argout] = fh(strucinput.(char(OPnames(i))),argout{:});
    end
    if Switch == 1
        strucoutput.(char(OPnames(i))) = strucinput.(char(OPnames(i)));
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

