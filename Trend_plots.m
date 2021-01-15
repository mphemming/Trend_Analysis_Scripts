%% Load in data

addpath(genpath('C:\Users\mphem\Documents\Work\UNSW\'))

options.data_dir = 'C:\Users\mphem\Documents\Work\UNSW\Trends\Data\';
options.plot_dir = 'C:\Users\mphem\Documents\Work\UNSW\Trends\Plots\';

% NRSPHB
NRSPHB_data = load([options.data_dir,'NRSPHB_data']);
NRSPHB_trends = load([options.data_dir,'NRSPHB_trends']);
%CH100



% Sort out time
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
tr = NRSPHB_trends.EEMD_trend{1};
p1 = plot(NRSPHB_trends.EEMD_t_conv(1).t,NRSPHB_trends.EEMD_T{1}-tr(1),'LineWidth',2)
p2 = plot(NRSPHB_trends.EEMD_t_conv(1).t,tr-tr(1),'LineWidth',2,'Color','r')
check_sig = tr-tr(1) > 0.4;
p3 = plot(NRSPHB_trends.EEMD_t_conv(1).t(check_sig),tr(check_sig)-tr(1),'LineWidth',2,'Color','k')

[p3a, p3b] = boundedline(NRSPHB_trends.EEMD_t_conv(1).t,zeros(size(NRSPHB_trends.EEMD_t_conv(1).t)),0.4)
set(p3b,'FaceColor','r','FaceAlpha',0.2)
set(p3a,'LineStyle','None')

set(gca,'LineWidth',2,'Box','On','FontSize',18,'YLim',[-3 4],'XLim',[datenum(1950,01,01) datenum(2021,01,01)],...
    'XTick',[datenum(1950,01,01) datenum(1960,01,01) datenum(1970,01,01) datenum(1980,01,01) datenum(1990,01,01) ...
    datenum(2000,01,01) datenum(2010,01,01) datenum(2020,01,01)],'XTickLabels',[{'1950'} {'1960'} {'1970'} {'1980'} {'1990'} {'2000'} {'2010'} {'2020'}])

ylabel('Temperature Anomaly [^\circC]');

leg = legend([p1 p2 p3 p3b],'De-seasoned temperatures','Non-significant EEMD Trend','Significant EEMD Trend','Boundary of non-significance');
set(leg,'Location','NorthWest','Box','Off','Position',[0.101466049382716 0.767275373893183 0.238859958518986 0.146990744076638]);

axes('Parent',gcf,...
    'Position',[0.620804398148148 0.11 0.284195601851853 0.815]);

p1 = plot(NRSPHB_trends.EEMD_t_conv(1).t,NRSPHB_trends.EEMD_T{1}-tr(1),'LineWidth',2)
hold on;

imf = NRSPHB_trends.EEMD_imfs{1};
a = 0;
for n = 1:size(imf,1)
    if n < 2
        a = a+4;
        p2 = plot(NRSPHB_trends.EEMD_t_conv(1).t,imf(n,:)-a,'LineWidth',2,'Color','k')
    else
        a = a+2;
        if n < size(imf,1)-1
            p2 = plot(NRSPHB_trends.EEMD_t_conv(1).t,imf(n,:)-a,'LineWidth',2,'Color','k')
        else
            p3 = plot(NRSPHB_trends.EEMD_t_conv(1).t,imf(n,:)-a,'LineWidth',2,'Color','r')
        end
    end
end

set(gca,'LineWidth',2,'Box','On','FontSize',18,'YLim',[-24 5],'XLim',[datenum(1950,01,01) datenum(2021,01,01)],...
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

leg = legend([p1 p2 p3],'De-seasoned temperatures','IMFs','IMFs used for trend');
set(leg,'Location','SouthWest','Box','Off');

print(gcf, '-dpng','-r400', [options.plot_dir,'EEMD_example_NRSPHB_2m'])

%% Example ITA plot

lin_1 = (0:0.026:0.5) + ones(1,20)*10;
lin_2 = (0:0.026:0.5)*-1  + ones(1,20)*10;
e_1 = (exp(0:0.052:1)-1)  + ones(1,20)*10;
e_2 = (exp(0:0.052:1)*-1) + ones(1,20)*10 + 1;

figure('units','normalized','position',[0 0.1 .9 .7]);

% Create axes
axes('Parent',gcf,...
    'Position',[0.13 0.11 0.382008101851852 0.814999999999993]);

p1 = plot(lin_1,'LineWidth',2);
hold on;
p1 = plot(lin_2,'LineWidth',2);
p1 = plot(e_1,'LineWidth',2);
p1 = plot(e_2,'LineWidth',2);

% Create axes
axes('Parent',gcf,...
    'Position',[0.528790509259259 0.11 0.376209490740741 0.81499999999999]);

hold on;
% lin_1
s1 = scatter(sort(lin_1(1:10)),sort(lin_1(11:20)),'filled')
s2 = scatter(sort(lin_2(1:10)),sort(lin_2(11:20)),'filled')
s3 = scatter(sort(e_1(1:10)),sort(e_1(11:20)),'filled')
s4 = scatter(sort(e_2(1:10)),sort(e_2(11:20)),'filled')
s5 = plot(9:11,9:11,'LineWidth',2,'Color','k');

%% ITA plots at CH100 and BMP120







