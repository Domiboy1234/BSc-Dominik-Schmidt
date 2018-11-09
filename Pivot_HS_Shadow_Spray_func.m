function [Headline,Export] = Pivot_HS_Shadow_Spray_func(strucinput)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
GET_fieldnames = fieldnames(strucinput);
Headline = [];

Doku_Nr = length(strucinput.(char(GET_fieldnames{1})).Docu_Names);

for i = 1 : Doku_Nr
    Headline = [Headline,{strucinput.(char(GET_fieldnames{1})).Docu_Names{i}}];
end
Headline = [Headline,{'shot_nr'},{'time'},{'total penetration'},{'spray direction'},{'spray angle 2 planes'},{'spray area'},{'split angle'},{'plumeangle_circle'}];

for i = 1 : length(GET_fieldnames)
    for j = 1 :size(strucinput.(char(GET_fieldnames{i})).total_penetration,2)
        time = strucinput.(char(GET_fieldnames{i})).time;
        Export_sub = cell(length(time),length(Headline));
        for k = 1 : Doku_Nr
            Export_sub(:,k)=strucinput.(char(GET_fieldnames{i})).Docu_Values(1,k);
        end
        Export_sub(:,Doku_Nr + 1) = {j};
        Export_sub(:,Doku_Nr + 2)=num2cell(time');
        Export_sub(:,Doku_Nr + 3)=num2cell(strucinput.(char(GET_fieldnames{i})).total_penetration(:,j));
        Export_sub(:,Doku_Nr + 4)=num2cell(strucinput.(char(GET_fieldnames{i})).spraydirection(:,j));
        Export_sub(:,Doku_Nr + 5)=num2cell(strucinput.(char(GET_fieldnames{i})).plumeangle_2planes(:,j));
        Export_sub(:,Doku_Nr + 6)=num2cell(strucinput.(char(GET_fieldnames{i})).sprayarea(:,j));
        Export_sub(:,Doku_Nr + 7)=num2cell(strucinput.(char(GET_fieldnames{i})).plumeangle_maxspraywidth(:,j));
        Export_sub(:,Doku_Nr + 8)=num2cell(strucinput.(char(GET_fieldnames{i})).plumeangle_circle(:,j));

%         try
%             Export_sub(:,16) = (strucinput.(char(GET_fieldnames{i})).GET_Parameters_Spray.Condition_Values(1,9));
%         catch
%              Export_sub(:,15) = {'210 kHz'};
%         end
        if i == 1 && j==1
            Export = Export_sub;
        else
            Export = [Export;Export_sub];
        end
    end
end
end

