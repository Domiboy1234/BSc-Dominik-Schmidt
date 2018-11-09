function [strucoutput] = ImportData_func(strucinput,current_path)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
EvaluationMode_Path = [pwd,'\Evaluation Mode'];
addpath(EvaluationMode_Path);

EvaluationMode_Path = dir(EvaluationMode_Path);
[indx,tf] = listdlg('PromptString','Select an Evaluation Mode:','SelectionMode','single','ListString',{EvaluationMode_Path.name});
fh = str2func(strrep(EvaluationMode_Path(indx).name,'.m',''));
OPnames = fieldnames(strucinput);
strucoutput = struct;
for i=1:length(OPnames)

    strucoutput.(char(OPnames{i})) = fh(strucinput.(char(OPnames{i})),'loaddata','Yes','current_path',current_path);

end
end

