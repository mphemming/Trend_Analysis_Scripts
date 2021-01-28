%% Load in data

addpath(genpath('C:\Users\mphem\Documents\Work\UNSW\'))
options.data_dir = 'C:\Users\mphem\Documents\Work\UNSW\Trends\Data\';
options.plot_dir = 'C:\Users\mphem\Documents\Work\UNSW\Trends\Plots\';

% NRSPHB
NRSPHB_data = load([options.data_dir,'NRSPHB_data']);
NRSPHB_trends = load([options.data_dir,'NRSPHB_trends']);
NRSPHB_trends_server = load([options.data_dir,'NRSPHB_trends_server']);
% NRSMAI
NRSMAI_data = load([options.data_dir,'NRSMAI_data']);
NRSMAI_trends = load([options.data_dir,'NRSMAI_trends']);
NRSMAI_trends_server = load([options.data_dir,'NRSMAI_trends_server']);
%CH100
CH100_data = load([options.data_dir,'CH100_data']);
CH100_trends = load([options.data_dir,'CH100_trends']);
CH100_trends_server = load([options.data_dir,'CH100_trends_server']);
%BMP120
BMP120_data = load([options.data_dir,'BMP120_data']);
BMP120_trends = load([options.data_dir,'BMP120_trends']);
BMP120_trends_server = load([options.data_dir,'BMP120_trends_server']);


%% Sort out time
% NRSPHB
% time
nt = size(NRSPHB_data.t,2);
for nn = 1:9
    for t = 1:nt
        a = squeeze(vertcat(NRSPHB_data.t(nn,t,:)))';
        NRSPHB_data.t_conv(nn).t(t) = datenum(convertCharsToStrings(a));
    end
end
% EEMD time
for nn = 1:9
    a = NRSPHB_trends.EEMD_t{nn};
    for t = 1:size(a,1)
        b = a(t,:);
        NRSPHB_trends.EEMD_t_conv(nn).t(t) = datenum(convertCharsToStrings(b));
    end
end

%% calculate average trends

%% Port Hacking

multiplier = 12*10;
[trends] = get_trends(NRSPHB_data, NRSPHB_trends,NRSPHB_trends_server,multiplier)



