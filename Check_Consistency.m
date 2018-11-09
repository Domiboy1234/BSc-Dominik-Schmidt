function error = Check_Consistency(smartstruct)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
fnames = fieldnames(smartstruct);
error = 0;
for i = 1 : length(fnames)
    NrofOP = length(smartstruct.(char(fnames{i})).folder);
    if ~(size(smartstruct.(char(fnames{i})).imagenames,1) == size(smartstruct.(char(fnames{i})).folder,1))
        disp(['Error in ',smartstruct.(char(fnames{i})).folder{1}])
        error = 1;
    end
end
end

