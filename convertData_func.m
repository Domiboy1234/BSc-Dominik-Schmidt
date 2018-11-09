function [strucoutput] = convertData_func(strucinput)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
OPnames = fieldnames(strucinput);
strucoutput = struct;
for i=1:length(OPnames)
    
    strucoutput.(char(OPnames{i})) = Export_Data(strucinput.(char(OPnames{i})));
end
end

