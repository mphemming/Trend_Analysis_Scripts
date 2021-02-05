%% Load in data

addpath(genpath('C:\Users\mphem\Documents\Work\UNSW\'))
options.data_dir = 'C:\Users\mphem\Documents\Work\UNSW\Trends\Data\';
options.plot_dir = 'C:\Users\mphem\Documents\Work\UNSW\Trends\Plots\';

% NRSPHB
NRSPHB_data = load([options.data_dir,'NRSPHB_data']);
NRSPHB_data_server = load([options.data_dir,'NRSPHB_data_server']);
NRSPHB_trends = load([options.data_dir,'NRSPHB_trends']);
NRSPHB_trends_server = load([options.data_dir,'NRSPHB_trends_server']);
NRSPHB_Salinity_analysis = load([options.data_dir,'NRSPHB_Salinity_analysis']);
% NRSMAI
NRSMAI_data = load([options.data_dir,'NRSMAI_data']);
NRSMAI_data_server = load([options.data_dir,'NRSMAI_data_server']);
NRSMAI_trends = load([options.data_dir,'NRSMAI_trends_server']); % using server for both for now because only 2 depths at moment
NRSMAI_trends_server = load([options.data_dir,'NRSMAI_trends_server']);
NRSMAI_Salinity_analysis = load([options.data_dir,'NRSMAI_Salinity_analysis']);

%% Sort out time

% NRSPHB

for nn = 1
    a = NRSPHB_trends.EEMD_t{nn};
    for t = 1:size(a,1)
        b = a(t,:);
        NRSPHB_trends.EEMD_t_conv(nn).t(t) = datenum(convertCharsToStrings(b));
    end
end

a = squeeze(NRSPHB_Salinity_analysis.EEMD_t_S(1,:,:));
for t = 1:size(a,1)
    b = a(t,:);
    NRSPHB_Salinity_analysis.EEMD_t_S_conv(nn).t(t) = datenum(convertCharsToStrings(b));
end

% NRSMAI

for nn = 1
    a = NRSMAI_trends.EEMD_t{nn};
    for t = 1:size(a,1)
        b = a(t,:);
        NRSMAI_trends.EEMD_t_conv(nn).t(t) = datenum(convertCharsToStrings(b));
    end
end

a = squeeze(NRSMAI_Salinity_analysis.EEMD_t_S(1,:,:));
for t = 1:size(a,1)
    b = a(t,:);
    NRSMAI_Salinity_analysis.EEMD_t_S_conv(nn).t(t) = datenum(convertCharsToStrings(b));
end

%% Create figure comparison

MAI_PSAL = NRSMAI_Salinity_analysis.EEMD_imfs_S.IMF_1(8,:);
MAI_tr = NRSMAI_trends_server.EEMD_imfs.IMF_1(8,:);
MAI_tr_t = NRSMAI_trends.EEMD_t_conv(1).t;

for n = 1:numel(NRSMAI_Salinity_analysis.EEMD_t_S_conv(1).t)
    check = MAI_tr_t >= NRSMAI_Salinity_analysis.EEMD_t_S_conv(1).t(n)-14 & ...
       NRSMAI_Salinity_analysis.EEMD_t_S_conv(1).t(n)+14;
   bin.t_MAI(n) = NRSMAI_Salinity_analysis.EEMD_t_S_conv(1).t(n);
   bin.PSAL_MAI(n) = MAI_PSAL(n);
   bin.TEMP_MAI(n) =nanmean(MAI_tr(check));
end
bin.PSAL_MAI = bin.PSAL_MAI-bin.PSAL_MAI(1);
bin.TEMP_MAI = bin.TEMP_MAI - bin.TEMP_MAI(1);

PHB_PSAL = NRSPHB_Salinity_analysis.EEMD_imfs_S.IMF_1(8,:);
PHB_tr = NRSPHB_trends_server.EEMD_imfs.IMF_1(8,:);
PHB_tr_t = NRSPHB_trends.EEMD_t_conv(1).t;

for n = 1:numel(NRSPHB_Salinity_analysis.EEMD_t_S_conv(1).t)
    check = PHB_tr_t >= NRSPHB_Salinity_analysis.EEMD_t_S_conv(1).t(n)-14 & ...
       NRSPHB_Salinity_analysis.EEMD_t_S_conv(1).t(n)+14;
   bin.t_PHB(n) = NRSPHB_Salinity_analysis.EEMD_t_S_conv(1).t(n);
   bin.PSAL_PHB(n) = PHB_PSAL(n);
   bin.TEMP_PHB(n) =nanmean(PHB_tr(check));
end
bin.PSAL_PHB = bin.PSAL_PHB-bin.PSAL_PHB(1);
bin.TEMP_PHB = bin.TEMP_PHB - bin.TEMP_PHB(1);
% Normalise
bin.PSAL_PHB = (bin.PSAL_PHB-nanmin(bin.PSAL_PHB))/(nanmax(bin.PSAL_PHB)-nanmin(bin.PSAL_PHB));
bin.PSAL_MAI = (bin.PSAL_MAI-nanmin(bin.PSAL_MAI))/(nanmax(bin.PSAL_MAI)-nanmin(bin.PSAL_MAI));
bin.TEMP_PHB = (bin.TEMP_PHB-nanmin(bin.TEMP_PHB))/(nanmax(bin.TEMP_PHB)-nanmin(bin.TEMP_PHB));
bin.TEMP_MAI = (bin.TEMP_MAI-nanmin(bin.TEMP_MAI))/(nanmax(bin.TEMP_MAI)-nanmin(bin.TEMP_MAI));

figure('units','normalized','position',[0 0.1 .6 .75]);

scatter(bin.PSAL_PHB,bin.TEMP_PHB);
hold on
scatter(bin.PSAL_MAI,bin.TEMP_MAI)

leg = legend('NRSPHB IMF_{A} 1993-2020','NRSMAI IMF_{A} 1993-2020');
set(leg,'Box','Off','Location','NorthWest');
% ylabel('NRSPHB Trend_{B+A} [^\circC]');
% ylabel('NRSMAI Trend_{B+A} [^\circC]');
ylim([0 1])
set(gca,'FontSize',18,'LineWidth',2,'Box','On');
grid on
ylabel('Normalised Temperature');
xlabel('Normalised Salinity');

% fit
c = isfinite(bin.TEMP_PHB) & isfinite(bin.PSAL_PHB);
[f_PHB,g_PHB,o_PHB] = fit(bin.PSAL_PHB(c)',bin.TEMP_PHB(c)','poly1');
c = isfinite(bin.TEMP_MAI) & isfinite(bin.PSAL_MAI);
[f_MAI,g_MAI,o_MAI] = fit(bin.PSAL_MAI(c)',bin.TEMP_MAI(c)','poly1');

print(gcf, '-dpng','-r400', [options.plot_dir,'PSAL_comparison'])

