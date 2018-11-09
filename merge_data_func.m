function [strucoutput] = merge_data_func(varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
running_number = 0;
strucoutput = struct;
for j = 1:length(varargin)
    strucinput = varargin{j};
    strucinput_fieldnames = fieldnames(strucinput);
    
    
    for i = 1 : length(strucinput_fieldnames)
        running_number = running_number +1;
        OPname = strucinput_fieldnames{i};
        OPname = strsplit(OPname,'_');
        newOPname = ['OP',num2str(running_number),'_',OPname{2}];
        strucoutput.(char(newOPname)) = strucinput.(char(strucinput_fieldnames{i}));
    end
    
    
    
    
end

