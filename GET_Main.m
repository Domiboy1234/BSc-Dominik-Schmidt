addpath([pwd '\VFunct\']); % Path for Functions of GET


Current_path = ('E:'); % Harddrive
Current_folder = ('Datensätze\Highspeed-Shadowgraphy'); % Folder on the Harddrive that has to be evaluated
Head_folder = ('Highspeed-Shadowgraphy'); % Part of the Current Folder, Main Folder of the Dataset
targetdata=struct;
targetdata(1).name='.tif';
emptystruct = struct;


%% Smartstructure finds all Data in the Current_folder
[ strucoutput_smart ] = smartstructure_func( Current_path,Current_folder, Head_folder,targetdata,emptystruct,1 );
% error = Check_Consistency(strucoutput_smart); % Consistency check if all single shots have the same number of images
% if error == 1
%     stop=1;
% end

%% Parameter Input for Evaluation
% template = struct;
% load([pwd,'\template']);
% [strucoutput_Para] = Parameter_func(strucoutput_smart,'oneforall',1,'template',template); % oneforall = 1 -> first set of parameters is used for the complete dataset, template -> no manual input
[strucoutput_Para] = Parameter_func(strucoutput_smart);
%% Evaluation
[strucoutput_Eval] = Evaluation_func(strucoutput_Para);

%% Export of Data to Excel and .mat - files
[Headline,Export] = Pivot_time_HS_Wall_Impingement_func(strucoutput_Eval);
mkdir([Current_path,'\',Head_folder,'Excel Export']);
xlswrite([Current_path,'\',Head_folder,'Excel Export','\Injektor_time.xlsx'],[Headline;Export],'Injektor_time','A1');

[Headline,Export] = Pivot_HS_Wall_Impingement_func(strucoutput_Eval);
mkdir([Current_path,'\',Head_folder,'Excel Export']);
xlswrite([Current_path,'\',Head_folder,'Excel Export','\Injektor.xlsx'],[Headline;Export],'Injektor','A1');

[strucoutput_Eval_Export] = ExportData_func(strucoutput_Eval);
SaveData_func(strucoutput_Eval_Export,Current_path,'Injektor');

%% for-loop for changing single properties of a structure containing Data
generic_struct = strucoutput_Eval;
fnames = fieldnames(generic_struct);
value = 'bwarea'; % name of the variable that has to be changed
for i = 1 : length(fnames)
generic_struct.(char(fnames{i})).(char(value)) = 1;
end

%% Usefull GET functions
strucinput = strucoutput_Eval;
[strucoutput] = Filter_func(strucinput); % Filter OPs for Documented Values, for example OP Number, Chamber pressure etc...
[strucoutput] = merge_data_func(strucinput1,strucinput2); % Merge multiple Datasets
[strucoutput] = LoadData_func([Current_path,'\Injektor10']); % Loads data saved by SaveData_func
[strucoutput] = ImportData_func(strucoutput,Current_path); % Creates a Class Object from loaded Data, current path being the HarddriveC
[strucoutput] = ExportData_func(strucinput);
i = 1;
[strucoutput] = struct(strucinput.(char(fnames{i}))); % makes a structure from a class element for examination of hidden properties


%% Usefull Wall_Impingement functions
make_Video(strucinput); % Creates Videos and Images with Edge-Overlay

%test
