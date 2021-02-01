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

%% Create figure comparison

figure('units','normalized','position',[0 0.1 .6 .75]);


MAI_EKE = NRSMAI_EKE.EEMD_imfs(7,:) + NRSMAI_EKE.EEMD_imfs(8,:);
MAI_tr = NRSMAI_trends_server.EEMD_trend_EAC{1};
MAI_tr_t = NRSMAI_trends.EEMD_t_conv(1).t;

for n = 1:numel(NRSMAI_EKE.EEMD_t_conv(1).t)
    check = MAI_tr_t >= NRSMAI_EKE.EEMD_t_conv(1).t(n)-14 & ...
       NRSMAI_EKE.EEMD_t_conv(1).t(n)+14;
   bin.t(n) = NRSMAI_EKE.EEMD_t_conv(1).t(n);
   bin.EKE_MAI(n) = MAI_EKE(n);
   bin.TEMP_MAI(n) =nanmean(MAI_tr(check));
end
bin.EKE_MAI = bin.EKE_MAI-bin.EKE_MAI(1);
bin.TEMP_MAI = bin.TEMP_MAI - bin.TEMP_MAI(1);

PHB_EKE = NRSPHB_EKE.EEMD_imfs(7,:) + NRSPHB_EKE.EEMD_imfs(8,:);
PHB_tr = NRSPHB_trends_server.EEMD_trend_EAC{1};
PHB_tr_t = NRSPHB_trends.EEMD_t_conv(1).t;

for n = 1:numel(NRSPHB_EKE.EEMD_t_conv(1).t)
    check = PHB_tr_t >= NRSPHB_EKE.EEMD_t_conv(1).t(n)-14 & ...
       NRSPHB_EKE.EEMD_t_conv(1).t(n)+14;
   bin.t(n) = NRSPHB_EKE.EEMD_t_conv(1).t(n);
   bin.EKE_PHB(n) = PHB_EKE(n);
   bin.TEMP_PHB(n) =nanmean(PHB_tr(check));
end
bin.EKE_PHB = bin.EKE_PHB-bin.EKE_PHB(1);
bin.TEMP_PHB = bin.TEMP_PHB - bin.TEMP_PHB(1);

scatter(bin.EKE_PHB,bin.TEMP_PHB);
hold on
scatter(bin.EKE_MAI,bin.TEMP_MAI)






