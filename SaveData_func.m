function SaveData_func(strucinput,savestr,Name)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
 inputinfo = whos('strucinput');
 if inputinfo.bytes < 2*10^9
     save([savestr,'\',Name],'strucinput')
 else
     NrofParts = ceil(inputinfo.bytes/(2*10^9));
     IndexParts = struct;
     fnames = fieldnames(strucinput);
     separatestruct = struct;
     sizepart = 0;
     OP=0;
     prtNr = 0;
     while length(fnames)>=OP
         prtNr=prtNr+1;
         sizepart = 0;
         while sizepart < 2*10^9
             OP = OP+1;
             if length(fnames)>=OP
             separatestruct.(char(['prt',num2str(prtNr)])).(char(fnames{OP})) = strucinput.(char(fnames{OP}));
             prt = separatestruct.(char(['prt',num2str(prtNr)]));
             prtinfo= whos('prt');
             sizepart = prtinfo.bytes;
             else
                 break
             end
         end
         
     end
     prtnames = fieldnames(separatestruct);
     Indexcell = cell(size(prtnames));
     for i = 1 : length(prtnames)
         prt = separatestruct.(char(prtnames{i}));
         save([savestr,'\',Name,'_prt',num2str(i)],'prt');
         Indexcell{i} = [Name,'_prt',num2str(i)];
     end
     IndexParts.Index = Indexcell;
     IndexParts.folder = savestr;
     save([savestr,'\',Name],'IndexParts');
 end
end

