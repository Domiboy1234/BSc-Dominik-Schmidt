function [ typen ] = dataname_func( typen,removechars )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
firstreplace = {'1';'2';'3';'4';'5';'6';'7';'8';'9';'0'};

for j = 1 : size(typen,1)
    temp = char(typen(j,1));
    %             temp2 = cell(size(temp,2),1);
    k = size(temp,2);
    l = 1;
    while k > 0
        u = strcmp(temp(1,l),removechars);
        if sum(u)==0
            l=l+1;
            k=k-1;
        else
            temp(l)='';
            k=k-1;
        end
    end
    
    typen(j,1)={temp};
end

for j = 1 : size(typen,1)
    if sum(strcmp(typen{j}(1),firstreplace))>0
        typen{j} = ['A' typen{j}];
    end
end

end

