function strucoutput = LoadData_func(loadstr)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
target_struct = load(loadstr);
if ~isfield(target_struct,'IndexParts')
     temp = load(loadstr);
     fnames = fieldnames(temp);
     strucoutput = temp.(char(fnames{1}));
else
    IndexParts = target_struct.IndexParts;
    strucoutput = struct;
    
    for i = 1 : size(IndexParts.Index,1)
        newpath = strsplit(loadstr,'\');
        newpath(length(newpath)) = IndexParts.Index(i);
        newpath = strjoin(newpath,'\');
        load(newpath);
        fnames = fieldnames(prt);
        for j = 1 : length(fnames)
            strucoutput.(char(fnames{j})) = prt.(char(fnames{j}));
        end
    end
end

