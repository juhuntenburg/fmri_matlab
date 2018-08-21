clear all
close all

%% Initialize variables.
filename = '180704C.txt';
cycles = 5;
nostim = 40;
stim = 20;
delimiter = ',';%'\t';
% if nargin<=3
    startRow    = 5;
    endRow      = inf;
% end

%% Format string for each line of text:
% be careful with the format, edit here for double (%f) or cell array of strings (%s), depending on the 
% data format in the columns of interest. For more information, see the TEXTSCAN documentation.

% formatSpec = '%f%n%n%f%f%f%f%*s%*s%*s%*s%*s%*s%[^\n\r]';
formatSpec = '%f%n%n%s%f%s%f%*s%*s%*s%*s%*s%*s%[^\n\r]';
% to import all columns (double format)
% formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]'; 

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string. 
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Allocate imported array to column variable names
time = dataArray{:, 1};
resp_rate = dataArray{:, 2};
resp_period = dataArray{:, 3};
temp = dataArray{:, 4};
event = dataArray{:, 5};
SpO2 = dataArray{:, 6};
temp_new = dataArray{:,7};

%% save in new mat file
monit.time = time; 
monit.resp_rate = resp_rate;
monit.resp_period = resp_period;
monit.temp = temp;
monit.event = event;
monit.SpO2 = SpO2;
monit.temp_new = temp_new;

%SAVE_DIR = '/Users/MFnew/Dropbox (Mainen Lab)/Mouse data_Madalena (1)/Ratbase_MouseLaser/';
%save([SAVE_DIR mouse filesep mouse '_' date '_monit.mat'],'monit');
DIR = [pwd filesep];
[f file ~] = fileparts(filename);
save([file '.mat'],'monit');

%%
plot(monit.event)

if length(find(monit.event>0))>1
    ttl_step_up = find(monit.event>0);
    k=1;
    for j = 2:cycles:length(ttl_step_up)
        difference=round(mean(diff(ttl_step_up(j:j+cycles-1))));
        r_rate{k} = resp_rate(ttl_step_up(j)-round(difference*nostim/(nostim+stim)):ttl_step_up(j+cycles-1)+difference-1);
        tem{k} = temp(ttl_step_up(j)-round(difference*nostim/(nostim+stim)):ttl_step_up(j+cycles-1)+difference-1);
        k=k+1;
    end
    save([file 'r_rate.mat'],'r_rate');
    save([file 'tem.mat'],'tem');
    k=1;
    if ~isempty(find(monit.SpO2>0,1))
        for j = 2:cycles:length(ttl_step_up)
            difference=round(mean(diff(ttl_step_up(j:j+cycles-1))));
            SpO{k} = SpO2(ttl_step_up(j)-round(difference*nostim/(nostim+stim)):ttl_step_up(j+cycles-1)+difference-1);
            k=k+1;
        end
        save([file 'SpO.mat'],'SpO');
    end 
end
%end