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
NRSPHB_SYDAP_analysis = load([options.data_dir,'SYDAP_TEMP_analysis']);
NRSPHB_EKE = load([options.data_dir,'NRSPHB_EKE_analysis']);

%% Sort out time

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

a = squeeze(NRSPHB_SYDAP_analysis.EEMD_t(1,:,:));
for t = 1:size(a,1)
    b = a(t,:);
    NRSPHB_SYDAP_analysis.EEMD_t_conv(nn).t(t) = datenum(convertCharsToStrings(b));
end

a = squeeze(NRSPHB_EKE.EEMD_t(1,:,:));
for t = 1:size(a,1)
    b = a(t,:);
    NRSPHB_EKE.EEMD_t_conv(nn).t(t) = datenum(convertCharsToStrings(b));
end


%% Figure to compare trends


% trend EAC

SYDAP_tr_2 = NRSPHB_SYDAP_analysis.EEMD_imfs(8,:)+NRSPHB_SYDAP_analysis.EEMD_imfs(7,:)+NRSPHB_SYDAP_analysis.EEMD_imfs(6,:);
SYDAP_tr = NRSPHB_SYDAP_analysis.EEMD_imfs(8,:)+NRSPHB_SYDAP_analysis.EEMD_imfs(7,:);
SYDAP_tr_t = NRSPHB_SYDAP_analysis.EEMD_t_conv(1).t;
PHB_PSAL = NRSPHB_Salinity_analysis.EEMD_imfs_S.IMF_1(8,:)+NRSPHB_Salinity_analysis.EEMD_imfs_S.IMF_1(9,:);
PSAL_tr_t = NRSPHB_Salinity_analysis.EEMD_t_S_conv(1).t;
PHB_tr = NRSPHB_trends_server.EEMD_imfs.IMF_1(8,:)+NRSPHB_trends_server.EEMD_imfs.IMF_1(7,:);
PHB_tr_t = NRSPHB_trends.EEMD_t_conv(1).t;
PHB_EKE = NRSPHB_EKE.EEMD_imfs(7,:)+NRSPHB_EKE.EEMD_imfs(8,:);
PHB_EKE_tr_t = NRSPHB_EKE.EEMD_t_conv(1).t;

t_grid = datenum(1953,01,01):1:datenum(2020,01,01);
for n = 1:numel(t_grid)
    check_T = PHB_tr_t >= t_grid(n)-14 & ...
       PHB_tr_t < t_grid(n)+14;    
   check_S = PSAL_tr_t >= t_grid(n)-14 & ...
       PSAL_tr_t < t_grid(n)+14;
    check_SYDAP = SYDAP_tr_t >= t_grid(n)-14 & ...
       SYDAP_tr_t < t_grid(n)+14;   
    check_EKE = PHB_EKE_tr_t >= t_grid(n)-14 & ...
       PHB_EKE_tr_t < t_grid(n)+14;      
   bin.t_PHB(n) = t_grid(n);
   bin.PSAL_PHB(n) = nanmean(PHB_PSAL(check_S));
   bin.TEMP_PHB(n) =nanmean(PHB_tr(check_T));
   bin.SYDAP(n) = nanmean(SYDAP_tr(check_SYDAP));
   bin.SYDAP_2(n) = nanmean(SYDAP_tr_2(check_SYDAP));
   bin.EKE(n) = nanmean(PHB_EKE(check_EKE));   
end
bin.SYDAP = bin.SYDAP-bin.SYDAP(6211);
bin.SYDAP_2 = bin.SYDAP_2-bin.SYDAP_2(6211);
bin.TEMP_PHB = bin.TEMP_PHB-bin.TEMP_PHB(153);
% bin.PSAL_PHB = bin.PSAL_PHB-bin.PSAL_PHB(122);
% bin.EKE = bin.EKE-bin.EKE(14612);

% trend

SYDAP_tr_2 = NRSPHB_SYDAP_analysis.EEMD_imfs(8,:)+NRSPHB_SYDAP_analysis.EEMD_imfs(7,:);
SYDAP_tr = NRSPHB_SYDAP_analysis.EEMD_imfs(8,:);
SYDAP_tr_t = NRSPHB_SYDAP_analysis.EEMD_t_conv(1).t;
PHB_PSAL = NRSPHB_Salinity_analysis.EEMD_imfs_S.IMF_1(9,:);
PSAL_tr_t = NRSPHB_Salinity_analysis.EEMD_t_S_conv(1).t;
PHB_tr = NRSPHB_trends_server.EEMD_imfs.IMF_1(8,:);
PHB_tr_t = NRSPHB_trends.EEMD_t_conv(1).t;
PHB_EKE = NRSPHB_EKE.EEMD_imfs(8,:);
PHB_EKE_tr_t = NRSPHB_EKE.EEMD_t_conv(1).t;

t_grid = datenum(1953,01,01):1:datenum(2020,01,01);
for n = 1:numel(t_grid)
    check_T = PHB_tr_t >= t_grid(n)-14 & ...
       PHB_tr_t < t_grid(n)+14;    
   check_S = PSAL_tr_t >= t_grid(n)-14 & ...
       PSAL_tr_t < t_grid(n)+14;
    check_SYDAP = SYDAP_tr_t >= t_grid(n)-14 & ...
       SYDAP_tr_t < t_grid(n)+14;   
    check_EKE = PHB_EKE_tr_t >= t_grid(n)-14 & ...
       PHB_EKE_tr_t < t_grid(n)+14;      
   bin_mono.t_PHB(n) = t_grid(n);
   bin_mono.PSAL_PHB(n) = nanmean(PHB_PSAL(check_S));
   bin_mono.TEMP_PHB(n) =nanmean(PHB_tr(check_T));
   bin_mono.SYDAP(n) = nanmean(SYDAP_tr(check_SYDAP));
   bin_mono.SYDAP_2(n) = nanmean(SYDAP_tr_2(check_SYDAP));
   bin_mono.EKE(n) = nanmean(PHB_EKE(check_EKE));   
end
bin_mono.SYDAP = bin_mono.SYDAP-bin_mono.SYDAP(6211);
bin_mono.SYDAP_2 = bin_mono.SYDAP_2-bin_mono.SYDAP_2(6211);
bin_mono.TEMP_PHB = bin_mono.TEMP_PHB-bin_mono.TEMP_PHB(153);
% bin_mono.PSAL_PHB = bin_mono.PSAL_PHB-bin_mono.PSAL_PHB(122);
% bin_mono.EKE = bin_mono.EKE-bin_mono.EKE(14612);



%% create figure

figure('units','normalized','position',[0 0 .5 .9]);

colors = cbrewer('qual','Set2',4);
set(gcf,'Color','w');

% Create axes
axes('Parent',gcf,...
    'Position',[0.10703125 0.543981481481482 0.690625 0.386831275720165]);

plot(bin.t_PHB,bin.TEMP_PHB,'LineWidth',3,'Color',colors(1,:));
hold on;
SYD = inpaint_nans(bin.SYDAP);
SYD(1:nanmin(find(isfinite(bin.SYDAP)))) = NaN;
scatter(bin.t_PHB(1:700:end),SYD(1:700:end),20,'MarkerFaceColor',colors(1,:),'MarkerEdgeColor',colors(1,:),'Marker','d');
SYD = inpaint_nans(bin.SYDAP_2);
SYD(1:nanmin(find(isfinite(bin.SYDAP_2)))) = NaN;
scatter(bin.t_PHB(1:700:end),SYD(1:700:end),20,'MarkerFaceColor',colors(1,:),'MarkerEdgeColor',colors(1,:),'Marker','*');
xlim([datenum(1953,01,01) datenum(2020,01,01)]);
set(gca,'LineWidth',2,'Box','Off','YTick','','XTickLabels','','YColor','w','XTick', ...
    [datenum(1960,01,01) datenum(1970,01,01) datenum(1980,01,01) datenum(1990,01,01) datenum(2000,01,01) datenum(2010,01,01)],...
    'XTickLabels','');
% datetick('x','YYYY','KeepLimits')
ylim([-0.3 1.6])
leg = legend('NRSPHB 2 m','Sydney Airport','Sydney Airport with additional IMF');
set(leg,'Location','NorthWest','Box','Off','Position',[0.199305555555556 0.845078873699748 0.303906255662441 0.0681584376857114]);


% Create axes
axes('Parent',gcf,...
    'Position',[0.0888020833333303 0.54295267489712 0.005 0.386831275720165]);
ylim([-0.3 1.6])
set(gca,'LineWidth',2,'FontSize',16,'YColor',colors(1,:),'XColor',colors(1,:));
ylabel('Temperature [^\circC]');


% Create axes
axes('Parent',gcf,...
    'Position',[0.10703125 0.543981481481482 0.690625 0.386831275720165]);

plot(bin.t_PHB,bin.PSAL_PHB,'LineWidth',3,'Color',colors(2,:));
ylim([-4E-3 8E-3])
set(gca,'Visible','Off')
xlim([datenum(1953,01,01) datenum(2020,01,01)]);

% Create axes
axes('Parent',gcf,...
    'Position',[0.855989583333331 0.544495884773663 0.005 0.386831275720165]);
ylim([-4E-3 8E-3])
set(gca,'LineWidth',2,'FontSize',14,'YColor',colors(2,:),'XColor',colors(2,:));
ylabel('Salinity');

% Create axes
axes('Parent',gcf,...
    'Position',[0.10703125 0.543981481481482 0.690625 0.386831275720165]);

plot(bin.t_PHB,bin.EKE,'LineWidth',3,'Color',colors(3,:));
ylim([0.006 0.018])
set(gca,'Visible','Off')
xlim([datenum(1953,01,01) datenum(2020,01,01)]);

% Create axes
axes('Parent',gcf,...
    'Position',[0.98 0.545010288065844 0.005 0.386831275720165]);

set(gca,'LineWidth',2,'FontSize',14,'YColor',colors(3,:),'XColor',colors(3,:),'YTick',[0.006 0.008 0.01 0.012 0.014 0.016 0.018],'YTickLabels',...
    [{'6'} {'8'} {'1'} {'1.2'} {'1.4'} {'1.6'} {'1.8'}]);
ylabel('EKE [ x 10^{-3} m^2 s^{-2}]','FontSize',14);
ylim([0.006 0.018])

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


% Create axes
axes('Parent',gcf,...
    'Position',[0.10703125 0.131430041152263 0.692708333333333 0.378600823045268]);

plot(bin_mono.t_PHB,bin_mono.TEMP_PHB,'LineWidth',3,'Color',colors(1,:));
hold on;
SYD = inpaint_nans(bin_mono.SYDAP);
SYD(1:nanmin(find(isfinite(bin_mono.SYDAP)))) = NaN;
scatter(bin_mono.t_PHB(1:700:end),SYD(1:700:end),20,'MarkerFaceColor',colors(1,:),'MarkerEdgeColor',colors(1,:),'Marker','d');
SYD = inpaint_nans(bin_mono.SYDAP_2);
SYD(1:nanmin(find(isfinite(bin_mono.SYDAP_2)))) = NaN;
scatter(bin.t_PHB(1:700:end),SYD(1:700:end),20,'MarkerFaceColor',colors(1,:),'MarkerEdgeColor',colors(1,:),'Marker','*');
xlim([datenum(1953,01,01) datenum(2020,01,01)]);
set(gca,'LineWidth',2,'Box','Off','YTick','','YColor','w', 'FontSize', 16, 'XTick', ...
    [datenum(1960,01,01) datenum(1970,01,01) datenum(1980,01,01) datenum(1990,01,01) datenum(2000,01,01) datenum(2010,01,01)],...
    'XTickLabels',[{'1960'} {'1970'} {'1980'} {'1990'} {'2000'} {'2010'}]);
ylim([-0.3 1.2])
leg = legend('NRSPHB 2 m','Sydney Airport','Sydney Airport with additional IMF');
set(leg,'Location','NorthWest','Box','Off','Position',[0.193055555555556 0.426354593790771 0.32994792294999 0.0712448575123838], ...
    'FontSize',10);

xlabel('Year');

% Create axes
axes('Parent',gcf,...
    'Position',[0.0893229166666634 0.131430041152263 0.00520833333333663 0.381172839506175]);
ylim([-0.3 1.2])
set(gca,'LineWidth',2,'FontSize',16,'YColor',colors(1,:),'XColor',colors(1,:));
ylabel('Temperature [^\circC]');


% Create axes
axes('Parent',gcf,...
    'Position',[0.10703125 0.131430041152263 0.692708333333333 0.378600823045268]);

plot(bin_mono.t_PHB,bin_mono.PSAL_PHB,'LineWidth',3,'Color',colors(2,:));
ylim([-2.75E-4 0.0001])
set(gca,'Visible','Off')
xlim([datenum(1953,01,01) datenum(2020,01,01)]);

% Create axes
axes('Parent',gcf,...
    'Position',[0.856510416666664 0.130401234567902 0.005 0.386831275720165]);

ylim([-2.75E-4 0.0001])
set(gca,'LineWidth',2,'FontSize',14,'YColor',colors(2,:),'XColor',colors(2,:));
ylabel('Salinity');

% Create axes
axes('Parent',gcf,...
    'Position',[0.10703125 0.131430041152263 0.692708333333333 0.378600823045268]);

plot(bin_mono.t_PHB,bin_mono.EKE,'LineWidth',3,'Color',colors(3,:));
ylim([1E-4 10E-4])
set(gca,'Visible','Off')
xlim([datenum(1953,01,01) datenum(2020,01,01)]);

% Create axes
axes('Parent',gcf,...
    'Position',[0.98 0.130401234567902 0.005 0.386831275720165]);

set(gca,'LineWidth',2,'FontSize',14,'YColor',colors(3,:),'XColor',colors(3,:),'YTick',[0.0002 0.0004 0.0006 0.0008 0.0010],'YTickLabels',...
    [{'2'} {'4'} {'6'} {'8'} {'10'}]);
ylabel('EKE [ x 10^{-4} m^2 s^{-2}]','FontSize',14);
ylim([1E-4 10E-4])



% Create textbox
annotation(gcf,'textbox',...
    [0.0978749999999999 0.878600823045268 0.0792083333333329 0.0421810699588483],...
    'String',{'(a)'},...
    'LineStyle','none',...
    'FontSize',24,...
    'FitBoxToText','off');

% Create textbox
annotation(gcf,'textbox',...
    [0.101520833333333 0.462448559670782 0.0792083333333326 0.0421810699588482],...
    'String','(b)',...
    'LineStyle','none',...
    'FontSize',24,...
    'FitBoxToText','off');

print(gcf, '-dpng','-r400', [options.plot_dir,'Trend_Sydney_comparison'])

