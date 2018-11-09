function [Headline,Export] = Pivot_time_HS_Wall_Impingement_func(strucinput)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
GET_fieldnames = fieldnames(strucinput);
Headline = [];

Doku_Nr = length(strucinput.(char(GET_fieldnames{1})).Docu_Names);

for i = 1 : Doku_Nr
    Headline = [Headline,{strucinput.(char(GET_fieldnames{1})).Docu_Names{i}}];
end
Headline = [Headline,{'shot_nr'},{'time'},{'bwarea_dark'},{'bwarea_bright'},{'bwarea_total'},{'time_step'},{'max_bwarea_dark'},{'max_bwarea_dark_warp'},{'max_time_dark'}];

for i = 1 : length(GET_fieldnames)
    for j = 1 :size(strucinput.(char(GET_fieldnames{i})).bwarea_dark,2)
        time = strucinput.(char(GET_fieldnames{i})).time;
        Export_sub = cell(length(time),length(Headline));
        for k = 1 : Doku_Nr
             Export_sub(:,k)={str2double(strucinput.(char(GET_fieldnames{i})).Docu_Values{1,k})};
        end
        Export_sub(:,Doku_Nr + 1) = {j};
        Export_sub(:,Doku_Nr + 2)=num2cell(time');
        Export_sub(:,Doku_Nr + 3)=num2cell(strucinput.(char(GET_fieldnames{i})).bwarea_dark(:,j));
        Export_sub(:,Doku_Nr + 4)=num2cell(strucinput.(char(GET_fieldnames{i})).bwarea_bright(:,j));
        Export_sub(:,Doku_Nr + 5)=num2cell(strucinput.(char(GET_fieldnames{i})).bwarea_total(:,j));
        Export_sub(:,Doku_Nr + 6)=num2cell(1/strucinput.(char(GET_fieldnames{i})).frequency);
        Export_sub(:,Doku_Nr + 7)=num2cell(strucinput.(char(GET_fieldnames{i})).max_bwarea_dark(1,j));
        Export_sub(:,Doku_Nr + 8)=num2cell(strucinput.(char(GET_fieldnames{i})).max_bwarea_dark_warp(1,j));
        Export_sub(:,Doku_Nr + 9)=num2cell(strucinput.(char(GET_fieldnames{i})).max_time_dark(1,j));

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

