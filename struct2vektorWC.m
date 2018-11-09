function [XY]=struct2vektorWC(str)
%x y aus einem Struct
[sp,anz]=size(str);
XY=zeros(sp,2);
for i=1:sp
    XY(i,:)=[str(i,1).WeightedCentroid];
end
