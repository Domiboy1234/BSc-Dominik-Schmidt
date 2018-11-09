function [strucoutput] = Evaluation_func(strucinput,varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% p = inputParser;
% p.KeepUnmatched = true;
% addRequired(p,'strucinput',@isstruct);
% addParameter(p,'ExportEdge',0);
% parse(p,strucinput,varargin{:});
% strucinput = p.Results.strucinput;
% ExportEdge = p.Results.ExportEdge;
EvaluationMode_Path = [pwd,'\Evaluation Mode'];
addpath(EvaluationMode_Path);
% EvaluationMode_Path = dir(EvaluationMode_Path);
% [indx,tf] = listdlg('PromptString','Select an Evaluation Mode:','SelectionMode','single','ListString',{EvaluationMode_Path.name});
% fh = str2func(strrep(EvaluationMode_Path(indx).name,'.m',''));
% Current_Data_Drive = uigetdir('C:\','Select current Data Path');
OPnames = fieldnames(strucinput);
% template_exists = 0;

for i=1:length(OPnames)
    [strucoutput.(char(OPnames(i)))] = Evaluation(strucinput.(char(OPnames(i))),varargin{:});
end
end

