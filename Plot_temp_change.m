%% Load in data

addpath(genpath('C:\Users\mphem\Documents\Work\UNSW\'))
options.data_dir = 'C:\Users\mphem\Documents\Work\UNSW\Trends\Data\';
options.plot_dir = 'C:\Users\mphem\Documents\Work\UNSW\Trends\Plots\';

% NRSPHB
NRSPHB_data = load([options.data_dir,'NRSPHB_data']);
NRSPHB_data_server = load([options.data_dir,'NRSPHB_data_server']);
NRSPHB_trends = load([options.data_dir,'NRSPHB_trends']);
NRSPHB_trends_server = load([options.data_dir,'NRSPHB_trends_server']);
% NRSMAI
NRSMAI_data = load([options.data_dir,'NRSMAI_data']);
NRSMAI_data_server = load([options.data_dir,'NRSMAI_data_server']);
NRSMAI_trends = load([options.data_dir,'NRSMAI_trends_server']); % using server for both for now because only 2 depths at moment
NRSMAI_trends_server = load([options.data_dir,'NRSMAI_trends_server']);
%CH100
CH100_data = load([options.data_dir,'CH100_data']);
CH100_trends = load([options.data_dir,'CH100_trends']);
CH100_trends_server = load([options.data_dir,'CH100_trends_server']);
%BMP120
BMP120_data = load([options.data_dir,'BMP120_data']);
BMP120_trends = load([options.data_dir,'BMP120_trends']);
BMP120_trends_server = load([options.data_dir,'BMP120_trends_server']);

MAI_depths = [2, 10, 20, 30, 40, 50, 85];
PHB_depths = [2, 19, 31, 40, 50, 59, 75, 81, 99];

%% Fixing depths where IMFs not capturing trend properly

a = NRSMAI_trends.EEMD_imfs.IMF_1;
NRSMAI_trends.EEMD_trend{1} = a(end,:) + a(end-1,:);
NRSMAI_trends.EEMD_trend_EAC{1} = a(end,:) + a(end-1,:) + a(end-2,:);

%% Get temp change

% NRSPHB
[Tchanges.NRSPHB.t1980_2020, Tchanges.NRSPHB.t1980_2020_EAC] = ...
    get_Tchange(NRSPHB_data, NRSPHB_data_server, NRSPHB_trends,NRSPHB_trends_server,1980,2020);
Tchanges.NRSPHB.t1980_2020_EAC([4,5,7,9]) = NaN;
% NRSMAI
[Tchanges.NRSMAI.t1980_2020, Tchanges.NRSMAI.t1980_2020_EAC] = ...
    get_Tchange(NRSMAI_data, NRSMAI_data_server, NRSMAI_trends,NRSMAI_trends_server,1980,2020);
Tchanges.NRSMAI.t1980_2020([4:5,7]) = NaN;
Tchanges.NRSMAI.t1980_2020_EAC([4,7]) = NaN;

%% Figure
figure('units','normalized','position',[0 0.1 .7 .8]);

% Create axes
axes('Parent',gcf,...
    'Position',[0.0846354166666667 0.11 0.405505952380952 0.815]);
x = 1:9;
plot(interp1(x(isfinite(Tchanges.NRSPHB.t1980_2020)),Tchanges.NRSPHB.t1980_2020(isfinite(Tchanges.NRSPHB.t1980_2020)), 1:9,'Linear'),...
    PHB_depths,'LineWidth',2)
hold on;
plot(interp1(x(isfinite(Tchanges.NRSPHB.t1980_2020_EAC)),Tchanges.NRSPHB.t1980_2020_EAC(isfinite(Tchanges.NRSPHB.t1980_2020_EAC)), 1:9,'Linear'),...
    PHB_depths,'LineWidth',2)
scatter(Tchanges.NRSPHB.t1980_2020,PHB_depths,'k','filled')
scatter(Tchanges.NRSPHB.t1980_2020_EAC,PHB_depths,'k','filled')

set(gca,'YDir','Reverse','LineWidth',2,'XGrid','Off','YGrid','Off','FontSize',16,'XLim',[0.2 1.2],'YLim',[0 85])
title('NRS Port Hacking');
xlabel('Temp. Change 1980-2020 [^\circC]')
ylabel('Depth [m]');

% Create axes
axes('Parent',gcf,...
    'Position',[0.507254464285714 0.11 0.397745535714286 0.815]);

x = 1:6;
p1 = plot(interp1(x(isfinite(Tchanges.NRSMAI.t1980_2020)),Tchanges.NRSMAI.t1980_2020(isfinite(Tchanges.NRSMAI.t1980_2020)), 1:6,'Linear'),...
    MAI_depths(1:end-1),'LineWidth',2)
hold on;
p2 = plot(interp1(x(isfinite(Tchanges.NRSMAI.t1980_2020_EAC)),Tchanges.NRSMAI.t1980_2020_EAC(isfinite(Tchanges.NRSMAI.t1980_2020_EAC)), 1:6,'Linear'),...
    MAI_depths(1:end-1),'LineWidth',2)
scatter(Tchanges.NRSMAI.t1980_2020,MAI_depths,'k','filled')
scatter(Tchanges.NRSMAI.t1980_2020_EAC,MAI_depths,'k','filled')

leg = legend([p2 p1],'Trend_{\rm{B+A}}','Trend_{\rm{B}}');
set(leg,'Location','SouthWest','Box','Off','FontSize',20);


set(gca,'YDir','Reverse','LineWidth',2,'YTickLabels','','XGrid','Off','YGrid','Off','FontSize',16,'XLim',[0.2 1.2],'YLim',[0 85],'XTick',[0.4 0.6 0.8 1 1.2])
title('NRS Maria Island');
xlabel('Temp. Change 1980-2020 [^\circC]')

% Create textbox
annotation(gcf,'textbox',...
    [0.0954940476190475 0.857638888888889 0.0525714285714285 0.0659722222222222],...
    'String',{'(a)'},...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',24,...
    'FitBoxToText','off');

% Create textbox
annotation(gcf,'textbox',...
    [0.518485119047619 0.855902777777778 0.0525714285714286 0.0659722222222222],...
    'String','(b)',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',24,...
    'FitBoxToText','off');

print(gcf, '-dpng','-r400', [options.plot_dir,'Tchange_1980_2020_NRS'])

