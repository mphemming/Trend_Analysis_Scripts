%% Load in data

addpath(genpath('C:\Users\mphem\Documents\Work\UNSW\'))
options.data_dir = 'C:\Users\mphem\Documents\Work\UNSW\Trends\Data\';
options.plot_dir = 'C:\Users\mphem\Documents\Work\UNSW\Trends\Plots\';

% NRSPHB
NRSPHB_data = load([options.data_dir,'NRSPHB_data']);
NRSPHB_data_server = load([options.data_dir,'NRSPHB_data_server']);
NRSPHB_trends = load([options.data_dir,'NRSPHB_trends']);
NRSPHB_trends_server = load([options.data_dir,'NRSPHB_trends_server']);
NRSPHB_EKE = load([options.data_dir,'NRSPHB_EKE_analysis']);
% NRSMAI
NRSMAI_data = load([options.data_dir,'NRSMAI_data']);
NRSMAI_data_server = load([options.data_dir,'NRSMAI_data_server']);
NRSMAI_trends = load([options.data_dir,'NRSMAI_trends_server']); % using server for both for now because only 2 depths at moment
NRSMAI_trends_server = load([options.data_dir,'NRSMAI_trends_server']);
NRSMAI_EKE = load([options.data_dir,'NRSMAI_EKE_analysis']);

% SOI

filename = 'C:\Users\mphem\Documents\Work\UNSW\Trends\Data\SOI.txt';
formatSpec = '%5s%7s%7s%7s%7s%7s%7s%7s%7s%7s%7s%7s%s%[^\n\r]';
fileID = fopen(filename,'r');

dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string',  'ReturnOnError', false);

fclose(fileID);
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2,3,4,5,6,7,8,9,10,11,12,13]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^[-/+]*\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end

R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

SOI = table;
SOI.VarName1 = cell2mat(raw(:, 1));
SOI.VarName2 = cell2mat(raw(:, 2));
SOI.VarName3 = cell2mat(raw(:, 3));
SOI.VarName4 = cell2mat(raw(:, 4));
SOI.VarName5 = cell2mat(raw(:, 5));
SOI.VarName6 = cell2mat(raw(:, 6));
SOI.VarName7 = cell2mat(raw(:, 7));
SOI.VarName8 = cell2mat(raw(:, 8));
SOI.VarName9 = cell2mat(raw(:, 9));
SOI.VarName10 = cell2mat(raw(:, 10));
SOI.VarName11 = cell2mat(raw(:, 11));
SOI.VarName12 = cell2mat(raw(:, 12));
SOI.VarName13 = cell2mat(raw(:, 13));

clearvars filename formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp R;


SOI_d.data(1,:) = [strct.VarName2];
SOI_d.data(2,:) = [strct.VarName3];
SOI_d.data(3,:) = [strct.VarName4];
SOI_d.data(4,:) = [strct.VarName5];
SOI_d.data(5,:) = [strct.VarName6];
SOI_d.data(6,:) = [strct.VarName7];
SOI_d.data(7,:) = [strct.VarName8];
SOI_d.data(8,:) = [strct.VarName9];
SOI_d.data(9,:) = [strct.VarName10];
SOI_d.data(10,:) = [strct.VarName11];
SOI_d.data(11,:) = [strct.VarName12];
SOI_d.data(12,:) = [strct.VarName13];

for n = 2:157
    SOI_d.Year(n) = SOI.VarName1(n);
    SOI_d.index(n) = nanmean(SOI_d.data(:,n))
end

SOI_d.Year = SOI_d.Year(76:end-2)
SOI_d.index = SOI_d.index(76:end-2)
SOI_d.time = datenum(SOI_d.Year,06,01);

%% Sort out time

% NRSPHB

for nn = 1
    a = NRSPHB_trends.EEMD_t{nn};
    for t = 1:size(a,1)
        b = a(t,:);
        NRSPHB_trends.EEMD_t_conv(nn).t(t) = datenum(convertCharsToStrings(b));
    end
end

a = squeeze(NRSPHB_EKE.EEMD_t(1,:,:));
for t = 1:size(a,1)
    b = a(t,:);
    NRSPHB_EKE.EEMD_t_conv(nn).t(t) = datenum(convertCharsToStrings(b));
end

% NRSMAI

for nn = 1
    a = NRSMAI_trends.EEMD_t{nn};
    for t = 1:size(a,1)
        b = a(t,:);
        NRSMAI_trends.EEMD_t_conv(nn).t(t) = datenum(convertCharsToStrings(b));
    end
end

a = squeeze(NRSMAI_EKE.EEMD_t(1,:,:));
for t = 1:size(a,1)
    b = a(t,:);
    NRSMAI_EKE.EEMD_t_conv(nn).t(t) = datenum(convertCharsToStrings(b));
end

%% bin data for comparison
a = NRSPHB_trends.EEMD_imfs.IMF_1;
A = NRSMAI_trends.EEMD_imfs.IMF_1;
for n = 1:numel(SOI_d.time)
    check = NRSPHB_trends.EEMD_t_conv(1).t > SOI_d.time(n)-datenum(0,6,01) &...
        NRSPHB_trends.EEMD_t_conv(1).t < SOI_d.time(n)+datenum(0,6,01);
    check_M = NRSMAI_trends.EEMD_t_conv(1).t > SOI_d.time(n)-datenum(0,6,01) &...
        NRSMAI_trends.EEMD_t_conv(1).t < SOI_d.time(n)+datenum(0,6,01);    
    aa = a(4,:)+a(5,:)+a(6,:);
    AA = A(4,:)+A(5,:)+A(6,:);
    bin.NRSPHB(n) = nanmean(aa(check));
    bin.NRSMAI(n) = nanmean(AA(check_M));
    bin.t(n) = SOI_d.time(n);
    bin.SOI(n) = SOI_d.index(n);
end

%% get fit
c = isfinite(bin.SOI) & isfinite(bin.NRSPHB);
[f,g,o] = fit(bin.SOI(c)',bin.NRSPHB(c)'*-1,'poly1');






