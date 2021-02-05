%% Load in data

addpath(genpath('C:\Users\mphem\Documents\Work\UNSW\'))
options.data_dir = 'C:\Users\mphem\Documents\Work\UNSW\Trends\Data\';
options.plot_dir = 'C:\Users\mphem\Documents\Work\UNSW\Trends\Plots\';

% NRSPHB
NRSPHB_data = load([options.data_dir,'NRSPHB_data']);
NRSPHB_data_server = load([options.data_dir,'NRSPHB_data_server']);
NRSPHB_trends = load([options.data_dir,'NRSPHB_trends']);
NRSPHB_trends_server = load([options.data_dir,'NRSPHB_trends_server']);

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

%% Example EEMD plot

figure('units','normalized','position',[0 0 .9 .9]);

axes('Parent',gcf,...
    'Position',[0.0646701388888889 0.11 0.52662037037037 0.815]);

hold on 
std_T = nanstd(NRSPHB_trends.EEMD_T{1});
tr = NRSPHB_trends.EEMD_trend{1} / std_T;
tr1 = tr(1)
tr = tr-tr1;
tr_EAC = NRSPHB_trends.EEMD_trend_EAC{1} / std_T;
tr_EAC = tr_EAC-tr_EAC(1);
p1 = plot(NRSPHB_trends.EEMD_t_conv(1).t,(NRSPHB_trends.EEMD_T{1}-tr1)/std_T,'LineWidth',2)
p2 = plot(NRSPHB_trends.EEMD_t_conv(1).t,tr,'LineWidth',2,'Color','r')
p3 = plot(NRSPHB_trends.EEMD_t_conv(1).t,tr_EAC,'LineWidth',2,'Color','r')
% determine where significant
for n = 1:numel(tr)
    t = NRSPHB_trends.EEMD_t_conv(1).t(n);
    f = find(NRSPHB_data.t_conv(1).t == t);
    if tr(n) <= NRSPHB_trends_server.EEMD_conf_std_limit(1,f)
        sig(n) = 0;
    else
        sig(n) = 1;
    end
    if tr_EAC(n) <= NRSPHB_trends_server.EEMD_conf_std_limit_EAC(1,f)
        sig_EAC(n) = 0;
    else
        sig_EAC(n) = 1;
    end     
end
date_sig = datestr(nanmin(NRSPHB_trends.EEMD_t_conv(1).t(sig == 1)));
p4 = plot(NRSPHB_trends.EEMD_t_conv(1).t(sig == 1),tr(sig == 1)-tr(1),'LineWidth',2,'Color','k')
p5 = plot(NRSPHB_trends.EEMD_t_conv(1).t(sig_EAC == 1),tr_EAC(sig_EAC == 1)-tr(1),'LineWidth',2,'Color',[1 .6 0])


[p3a, p3b] = boundedline(NRSPHB_data.t_conv(1).t,zeros(size(NRSPHB_data.t_conv(1).t)),NRSPHB_trends_server.EEMD_conf_std_limit_EAC(1,:))
[p4a, p4b] = boundedline(NRSPHB_data.t_conv(1).t,zeros(size(NRSPHB_data.t_conv(1).t)),NRSPHB_trends_server.EEMD_conf_std_limit(1,:))
set(p3b,'FaceColor','r','FaceAlpha',0.2); set(p3a,'LineStyle','None');
set(p4b,'FaceColor','r','FaceAlpha',0.3); set(p4a,'LineStyle','None');

set(gca,'LineWidth',2,'Box','On','FontSize',18,'YLim',[-3 4],'XLim',[datenum(1950,01,01) datenum(2021,01,01)],...
    'XTick',[datenum(1950,01,01) datenum(1960,01,01) datenum(1970,01,01) datenum(1980,01,01) datenum(1990,01,01) ...
    datenum(2000,01,01) datenum(2010,01,01) datenum(2020,01,01)],'XTickLabels',[{'1950'} {'1960'} {'1970'} {'1980'} {'1990'} {'2000'} {'2010'} {'2020'}])

ylabel('Temperature Anomaly [^\circC]');

leg = legend([p1 p4 p5 p3b p4b],'De-seasoned temperatures', ...
    'Significant EEMD Trend_{B}','Significant EEMD Trend_{B+A}','Boundary of insignificance Trend_{B+A}','Boundary of insignificance Trend_{B}');
set(leg,'Location','NorthWest','Box','Off','Position',[0.101466049382716 0.737275373893183 0.238859958518986 0.146990744076638],...
    'FontSize',14);

axes('Parent',gcf,...
    'Position',[0.620804398148148 0.11 0.284195601851853 0.815]);

p1 = plot(NRSPHB_trends.EEMD_t_conv(1).t,NRSPHB_trends.EEMD_T{1}-tr(1),'LineWidth',2)
hold on;

imf = NRSPHB_trends.EEMD_imfs.IMF_1;
a = 0;
for n = 1:size(imf,1)-1
    if n < 2
        a = a+4;
        p2 = plot(NRSPHB_trends.EEMD_t_conv(1).t,imf(n,:)-a,'LineWidth',2,'Color',[0 .4 .4])
    else
        a = a+2;
        if n < size(imf,1)-2
            p3 = plot(NRSPHB_trends.EEMD_t_conv(1).t,imf(n,:)-a,'LineWidth',2,'Color',[0 .4 .4])
        end
        if n == size(imf,1)-2
            p4 = plot(NRSPHB_trends.EEMD_t_conv(1).t,imf(n,:)-a,'LineWidth',2,'Color',[.6 .4 .8]) 
        end
        if n == size(imf,1)-1
            p5 = plot(NRSPHB_trends.EEMD_t_conv(1).t,imf(n,:)-a,'LineWidth',2,'Color','k')
        end
    end
end

p6 = plot(NRSPHB_trends.EEMD_t_conv(1).t,tr_EAC-tr(1)-20,'LineWidth',2,'Color',[1 .6 0])

set(gca,'LineWidth',2,'Box','On','FontSize',18,'YLim',[-28 5],'XLim',[datenum(1950,01,01) datenum(2021,01,01)],...
    'XTick',[datenum(1960,01,01) datenum(1980,01,01) datenum(2000,01,01) datenum(2020,01,01)],...
    'XTickLabels',[{'1960'} {'1980'} {'2000'} {'2020'}],'YTickLabels','','YTick','')

% Create textbox
annotation(gcf,'textbox',...
    [0.0669722222222224 0.870370370370371 0.0389305555555558 0.0421810699588474],...
    'String',{'(a)'},...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',18,...
    'FitBoxToText','off');

% Create textbox
annotation(gcf,'textbox',...
    [0.623974537037037 0.86985596707819 0.0389305555555558 0.0421810699588474],...
    'String','(b)',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',18,...
    'FitBoxToText','off');

leg = legend([p2 p4 p5 p6],'IMFs','IMF_{A}','Trend_{B}','Trend_{B+A}');
set(leg,'Location','SouthWest','Box','Off','FontSize',14);

print(gcf, '-dpng','-r400', [options.plot_dir,'EEMD_example_NRSPHB_2m'])



