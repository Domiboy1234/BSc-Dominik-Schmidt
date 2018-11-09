function [ strucoutput, layer ] = smartstructure_func( current_path, folder, HeadFolder, targetdata, strucinput, naminglayer )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

current_folder = dir([current_path,'\',folder]);
removechars = {'.';' ';'_';'-';'#';'[';']';'/';'\';'"';'{';'}';'(';')';'°';'ß';',';'='};
data_found =zeros(length(targetdata),1) ;
datavar=zeros(length(targetdata),1); %internal runnning variable for data aquisition
layer = NaN;



for i = 1:length(current_folder)
    %% Recursion for folder structure
    if current_folder(i).isdir==1 && ~strcmp(current_folder(i).name,'.')  && ~strcmp(current_folder(i).name,'..')
        strucname = dataname_func({current_folder(i).name},removechars);

        [ strucinputnew, layer ] = smartstructure_func(current_path, [folder '\' current_folder(i).name] ,HeadFolder,targetdata,strucinput, naminglayer );
        if layer == naminglayer
            fnames = fieldnames(strucinput);
            for k = 1:length(fnames)
                if isfield(strucinputnew,fnames{k})
                     strucinputnew = rmfield(strucinputnew,fnames{k});
                end
            end
            fnamesnew = fieldnames(strucinputnew);
            for k = 1:length(fnamesnew)
                strucinput.(char(['OP',num2str(length(fnames)+k),'_',fnamesnew{k}]))=strucinputnew.(char(fnamesnew{k}));
            end
        else
            strucinput = strucinputnew;
        end
        layer = layer +1;
        
        
    elseif current_folder(i).isdir==0
        %% Data Handling
        for j=1:length(targetdata)
            
            if ~isempty(strfind(current_folder(i).name,targetdata(j).name)) ==1
                layer = 0;
                if strcmp(targetdata(j).name,'.tif')==1
                    data_found(j,1) = 1;
                    datavar(j,1)=datavar(j,1)+1;
                    imagenames( datavar(j,1))={current_folder(i).name};
                elseif strcmp(targetdata(j).name,'.csv')==1
                    data_found(j,1) = 1;
                    data = energyread_func([current_path,'\',folder '\' current_folder(i).name]);
                    timecsv(:,datavar(j,1)) = data{:,1};
                    energycsv(:,datavar(j,1)) = data{:,2};
                else
                    %             strucinput.(char(strucname))=0;
                end
            end
        end
        
        
    end
    
end

for j=1:length(targetdata)
if  data_found(j,1) == 1 && strcmp(targetdata(j).name,'.tif')==1
   
    nameOPpre = strsplit(folder,'\');
    nameOPpre=nameOPpre(length(nameOPpre)-naminglayer);
    nameOPpre=dataname_func(nameOPpre,removechars);
    nameOP = nameOPpre{1};
%     nameOP=['OP',num2str(inputlength+1),'_',nameOPpre{1}];
    
    Pos_Nr = 1;
    if isfield(strucinput,nameOP)
%     while isfield(strucinput.(char(nameOP)),['POS',num2str(Pos_Nr),'_folder'])
%        Pos_Nr = Pos_Nr +1; 
%     end
        if isfield(strucinput.(char(nameOP)),'folder')
            Pos_Nr = 1 + length(strucinput.(char(nameOP)).folder);
        end
    end  
%     strucinput.(char(nameOP)).(char(['POS',num2str(Pos_Nr),'_folder']))(Pos_Nr) = {folder};
%     strucinput.(char(nameOP)).(char(['POS',num2str(Pos_Nr),'_imagenames']))(Pos_Nr,:) = imagenames;
    strucinput.(char(nameOP)).folder(Pos_Nr,1) = {folder};
    strucinput.(char(nameOP)).HeadFolder = HeadFolder;
    try
    strucinput.(char(nameOP)).imagenames(Pos_Nr,:) = imagenames;
    catch
        test=1
    end
     strucinput.(char(nameOP)).Pos_Nr = Pos_Nr;
     strucinput.(char(nameOP)).current_path(Pos_Nr,1) = {current_path};
elseif data_found(j,1) == 1 && strcmp(targetdata(j).name,'.csv')==1
    strucinput.(char(nameOP)).folder = folder;
    strucinput.(char(nameOP)).timecsv = timecsv;
    strucinput.(char(nameOP)).energycsv = energycsv;
    
end
end


strucoutput=strucinput;
end

