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

%% Remove some depths that lack data affecting trends

% NRSPHB
% missing 1950s/1960s with higher maxima than 1970s

NRSPHB_trends.EEMD_t{6} = [];
NRSPHB_trends.EEMD_T{6} = [];
NRSPHB_trends.EEMD_trend{6} = [];
NRSPHB_trends.EEMD_trend_EAC{6} = [];
NRSPHB_trends_server.EEMD_t{6} = [];
NRSPHB_trends_server.EEMD_T{6} = [];
NRSPHB_trends_server.EEMD_trend{6} = [];
NRSPHB_trends_server.EEMD_trend_EAC{6} = [];

NRSPHB_trends.EEMD_t{8} = [];
NRSPHB_trends.EEMD_T{8} = [];
NRSPHB_trends.EEMD_trend{8} = [];
NRSPHB_trends.EEMD_trend_EAC{8} = [];
NRSPHB_trends_server.EEMD_t{8} = [];
NRSPHB_trends_server.EEMD_T{8} = [];
NRSPHB_trends_server.EEMD_trend{8} = [];
NRSPHB_trends_server.EEMD_trend_EAC{8} = [];

%% Plot ITA for CH100 and BMP120

figure('units','normalized','position',[0 0.1 .8 .75]);

% Create axes
axes(gcf,'Position',[0.0818684895833333 0.11 0.397786458333333 0.814999999999983]);

cmap = cbrewer('qual','Paired',12);

clear s
for n = 1:11
  s(n) = scatter(CH100_trends_server.ITA_stats{n}.TEMP_half_1, CH100_trends_server.ITA_stats{n}.TEMP_half_2,10,'MarkerEdgeColor',cmap(n,:),'MarkerFaceColor',cmap(n,:))
  hold on
end
plot(-6:6,-6:6,'k','LineWidth',2);
xlim([-6 6]); ylim([-6 6]);
set(gca,'LineWidth',2,'Box','On','FontSize',16);
grid on
xlabel('Temperature Anomaly [^\circC] 2009 to 2015')
ylabel('Temperature Anomaly [^\circC] 2015 to 2021')
title('CH100');

clear ss
for n = 1:numel(s)
    if n == 1
        ss = s(n);
    else
        ss = [ss s(n)];
    end
end

leg = legend(ss,[{'10.5'}, {'20'}, {'27.5'}, {'35.5'}, {'43.5'}, {'51.5'}, {'59.5'}, {'67.5'}, {'75.5'}, {'84.5'}, {'91.5'}]);
set(leg,'Location','SouthEast','Box','Off');

% Create axes
axes(gcf,'Position',[0.526529947916667 0.11 0.404947916666667 0.814999999999977]);


for n = 1:12
  s(n) = scatter(BMP120_trends_server.ITA_stats{n}.TEMP_half_1, BMP120_trends_server.ITA_stats{n}.TEMP_half_2,10,'MarkerEdgeColor',cmap(n,:),'MarkerFaceColor',cmap(n,:))
  hold on
end
plot(-6:6,-6:6,'k','LineWidth',2);
xlim([-6 6]); ylim([-6 6]);
set(gca,'LineWidth',2,'Box','On','FontSize',16,'YTickLabels','');
grid on
xlabel('Temperature Anomaly [^\circC] 2011 to 2016')
ylabel('Temperature Anomaly [^\circC] 2016 to 2021')
title('BMP120');

for n = 1:numel(s)
    if n == 1
        ss = s(n);
    else
        ss = [ss s(n)];
    end
end

leg = legend(ss,[{'18.5'}, {'27.5'}, {'35'}, {'43'}, {'50.5'}, {'58.5'}, {'67'}, {'75'}, {'83.5'}, {'91.5'}, {'99.5'}, {'107.5'}]);
set(leg,'Location','SouthEast','Box','Off');

% Create textbox
annotation(gcf,'textbox',...
    [0.0921458333333333 0.84126984126984 0.0543385416666667 0.075396825396825],...
    'String',{'(c)'},...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',22,...
    'FitBoxToText','off');

% Create textbox
annotation(gcf,'textbox',...
    [0.535830729166667 0.841931216931211 0.0543385416666666 0.0753968253968249],...
    'String','(d)',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',22,...
    'FitBoxToText','off');

% Create textbox
annotation(gcf,'textbox',...
    [0.280296875 0.392857142857135 0.235653645833333 0.0753968253968245],...
    'Color',[0.650980392156863 0.650980392156863 0.650980392156863],...
    'String','Decreasing trend',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FitBoxToText','off');

% Create textbox
annotation(gcf,'textbox',...
    [0.727888020833334 0.39484126984126 0.235653645833333 0.0753968253968245],...
    'Color',[0.650980392156863 0.650980392156863 0.650980392156863],...
    'String','Decreasing trend',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FitBoxToText','off');

% Create textbox
annotation(gcf,'textbox',...
    [0.151065104166667 0.535052910052903 0.235653645833333 0.0753968253968249],...
    'Color',[0.650980392156863 0.650980392156863 0.650980392156863],...
    'String','Increasing trend',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FitBoxToText','off');

% Create textbox
annotation(gcf,'textbox',...
    [0.594098958333334 0.533068783068774 0.235653645833333 0.0753968253968249],...
    'Color',[0.650980392156863 0.650980392156863 0.650980392156863],...
    'String','Increasing trend',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FitBoxToText','off');

print(gcf, '-dpng','-r400', [options.plot_dir,'ITA_comparison_CH100_BMP120'])

%% Plot ITA for PHB and MAI

figure('units','normalized','position',[0 0.1 .8 .75]);

% Create axes
axes(gcf,'Position',[0.0818684895833333 0.11 0.397786458333333 0.814999999999983]);

cmap = cbrewer('qual','Paired',12);

clear s
for n = [1:5,7,9]
  s(n) = scatter(NRSPHB_trends_server.ITA_stats{n}.TEMP_half_1, NRSPHB_trends_server.ITA_stats{n}.TEMP_half_2,10,'MarkerEdgeColor',cmap(n,:),'MarkerFaceColor',cmap(n,:))
  hold on
end
plot(-6:6,-6:6,'k','LineWidth',2);
xlim([-4 4]); ylim([-4 4]);
set(gca,'LineWidth',2,'Box','On','FontSize',16);
grid on
xlabel('Temperature Anomaly [^\circC] 1953 to 1987')
ylabel('Temperature Anomaly [^\circC] 1987 to 2021')
title('NRSPHB');

clear ss
for n = [1:5,7,9]
    if n == 1
        ss = s(n);
    else
        ss = [ss s(n)];
    end
end

leg = legend(ss,[{'2'}, {'19'}, {'31'}, {'40'}, {'50'}, {'75'}, {'99'}]);
set(leg,'Location','SouthEast','Box','Off');

% Create axes
axes(gcf,'Position',[0.526529947916667 0.11 0.404947916666667 0.814999999999977]);

clear s
for n = [1,3,6]
  s(n) = scatter(NRSMAI_trends_server.ITA_stats{n}.TEMP_half_1, NRSMAI_trends_server.ITA_stats{n}.TEMP_half_2,10,'MarkerEdgeColor',cmap(n,:),'MarkerFaceColor',cmap(n,:))
  hold on
end
plot(-6:6,-6:6,'k','LineWidth',2);
xlim([-4 4]); ylim([-4 4]);
set(gca,'LineWidth',2,'Box','On','FontSize',16,'YTickLabels','');
grid on
xlabel('Temperature Anomaly [^\circC] 1945 to 1983')
ylabel('Temperature Anomaly [^\circC] 1983 to 2021')
title('NRSMAI');

clear ss
for n = [1,3,6]
    if n == 1
        ss = s(n);
    else
        ss = [ss s(n)];
    end
end

leg = legend(ss,[{'2'}, {'20'}, {'50'}, {'85'}]);
set(leg,'Location','SouthEast','Box','Off');

% Create textbox
annotation(gcf,'textbox',...
    [0.0921458333333333 0.84126984126984 0.0543385416666667 0.075396825396825],...
    'String',{'(a)'},...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',22,...
    'FitBoxToText','off');

% Create textbox
annotation(gcf,'textbox',...
    [0.535830729166667 0.841931216931211 0.0543385416666666 0.0753968253968249],...
    'String','(b)',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',22,...
    'FitBoxToText','off');

% Create textbox
annotation(gcf,'textbox',...
    [0.280296875 0.392857142857135 0.235653645833333 0.0753968253968245],...
    'Color',[0.650980392156863 0.650980392156863 0.650980392156863],...
    'String','Decreasing trend',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FitBoxToText','off');

% Create textbox
annotation(gcf,'textbox',...
    [0.727888020833334 0.39484126984126 0.235653645833333 0.0753968253968245],...
    'Color',[0.650980392156863 0.650980392156863 0.650980392156863],...
    'String','Decreasing trend',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FitBoxToText','off');

% Create textbox
annotation(gcf,'textbox',...
    [0.1 0.535052910052903 0.235653645833333 0.0753968253968249],...
    'Color',[0.650980392156863 0.650980392156863 0.650980392156863],...
    'String','Increasing trend',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FitBoxToText','off');

% Create textbox
annotation(gcf,'textbox',...
    [0.55 0.533068783068774 0.235653645833333 0.0753968253968249],...
    'Color',[0.650980392156863 0.650980392156863 0.650980392156863],...
    'String','Increasing trend',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FitBoxToText','off');

print(gcf, '-dpng','-r400', [options.plot_dir,'ITA_comparison_NRSPHB_NRSMAI'])
