function [Headline,Export] = PivotConti_HS_Shadow_Spray_func(strucinput)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
GET_fieldnames = fieldnames(strucinput);

Doku_Nr = length(strucinput.(char(GET_fieldnames{1})).Docu_Names);

Headline = [{'time step'},{'Index'},{'Slice (shot number)'},{'Time'},{'time2'},{'FUP'},{'Position'},{'Ti'},{'INJ'},{'REF'},{'Area'},{'Perim'},{'Mean'},...
            {'Sd'},{'Min'},{'Max'},{'XM'},{'YM'},{'JPen'},{'JetWidth'},{'spDir'},{'splitAngle'},{'splitAngleLeft'},{'splitAngleRight'},{'LeftAngle3'},{'LeftAngle10'},...
            {'LeftAngle30'},{'LeftAngle60'},{'RightAngle3'},{'RightAngle10'},{'RightAngle30'},{'RightAngle60'},{'spNearAngle'},{'p_chamber'},{'T_fuel'},{'Fuel'},{'Date'},{'Clocktime'}];
      
for i = 1 : length(GET_fieldnames)
    for j = 1 :size(strucinput.(char(GET_fieldnames{i})).total_penetration,2)
        time = strucinput.(char(GET_fieldnames{i})).time;
        Export_sub = cell(length(time),length(Headline));
      
        Export_sub(:,3) = {j};
        Export_sub(:,4) = num2cell(time');
        Export_sub(:,5) = num2cell(time');
        Export_sub(:,6) = strucinput.(char(GET_fieldnames{i})).Docu_Values(1,6);
        Export_sub(:,7) = strucinput.(char(GET_fieldnames{i})).Docu_Values(1,8);
        Export_sub(:,8) = {3000};
        Export_sub(:,9) = {[strucinput.(char(GET_fieldnames{i})).Docu_Values{1,1},strucinput.(char(GET_fieldnames{i})).Docu_Values{1,2},strucinput.(char(GET_fieldnames{i})).Docu_Values{1,3}]};
        Export_sub(:,11) = num2cell(strucinput.(char(GET_fieldnames{i})).sprayarea(:,j));
        Export_sub(:,19) = num2cell(strucinput.(char(GET_fieldnames{i})).total_penetration(:,j));
        Export_sub(:,21) = num2cell(strucinput.(char(GET_fieldnames{i})).spraydirection(:,j));
        Export_sub(:,22) = num2cell(strucinput.(char(GET_fieldnames{i})).plumeangle_maxspraywidth(:,j));
        Export_sub(:,33) = num2cell(strucinput.(char(GET_fieldnames{i})).plumeangle_2planes(:,j));
        Export_sub(:,34) =  strucinput.(char(GET_fieldnames{i})).Docu_Values(1,4);
        Export_sub(:,35) =  strucinput.(char(GET_fieldnames{i})).Docu_Values(1,7);
        Export_sub(:,36) =  strucinput.(char(GET_fieldnames{i})).Docu_Values(1,9);
        Export_sub(:,37) = {char(strucinput.(char(GET_fieldnames{i})).Date)};
        Export_sub(:,38) = {char(strucinput.(char(GET_fieldnames{i})).Time)};
        
%         Export_sub(:,Doku_Nr + 2)=num2cell(time');
%         Export_sub(:,Doku_Nr + 3)=num2cell(strucinput.(char(GET_fieldnames{i})).total_penetration(:,j));
%         Export_sub(:,Doku_Nr + 4)=num2cell(strucinput.(char(GET_fieldnames{i})).spraydirection(:,j));
%         Export_sub(:,Doku_Nr + 5)=num2cell(strucinput.(char(GET_fieldnames{i})).plumeangle_2planes(:,j));
%         Export_sub(:,Doku_Nr + 6)=num2cell(strucinput.(char(GET_fieldnames{i})).sprayarea(:,j));
%         Export_sub(:,Doku_Nr + 7)=num2cell(strucinput.(char(GET_fieldnames{i})).plumeangle_maxspraywidth(:,j));
%         Export_sub(:,Doku_Nr + 8)=num2cell(strucinput.(char(GET_fieldnames{i})).plumeangle_circle(:,j));

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

