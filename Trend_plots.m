%% Load in data

addpath(genpath('C:\Users\mphem\Documents\Work\UNSW\'))

options.data_dir = 'C:\Users\mphem\Documents\Work\UNSW\Trends\Data\';
options.plot_dir = 'C:\Users\mphem\Documents\Work\UNSW\Trends\Plots\';

% NRSPHB
NRSPHB_data = load([options.data_dir,'NRSPHB_data']);
NRSPHB_trends = load([options.data_dir,'NRSPHB_trends']);
NRSPHB_trends_server = load([options.data_dir,'NRSPHB_trends_server']);
% NRSMAI
% NRSMAI_data = load([options.data_dir,'NRSMAI_data']);
NRSMAI_trends = load([options.data_dir,'NRSMAI_trends']);
NRSMAI_trends_server = load([options.data_dir,'NRSMAI_trends_server']);
%CH100
CH100_trends = load([options.data_dir,'CH100_trends']);
%BMP120
BMP120_trends = load([options.data_dir,'BMP120_trends']);

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
std_T = nanstd(NRSPHB_trends.EEMD_T{1});
tr = NRSPHB_trends.EEMD_trend{1} / std_T;
tr = tr-tr(1);
p1 = plot(NRSPHB_trends.EEMD_t_conv(1).t,(NRSPHB_trends.EEMD_T{1}-tr(1))/std_T,'LineWidth',2)
p2 = plot(NRSPHB_trends.EEMD_t_conv(1).t,tr,'LineWidth',2,'Color','r')
% determine where significant
for n = 1:numel(tr)
    t = NRSPHB_trends.EEMD_t_conv(1).t(n);
    f = find(NRSPHB_data.t_conv(1).t == t);
    if tr(n) <= NRSPHB_trends_server.EEMD_conf_std_limit(1,f)
        sig(n) = 0;
    else
        sig(n) = 1;
    end
end
date_sig = datestr(nanmin(NRSPHB_trends.EEMD_t_conv(1).t(sig == 1)));
p3 = plot(NRSPHB_trends.EEMD_t_conv(1).t(sig == 1),tr(sig == 1)-tr(1),'LineWidth',2,'Color','k')

[p3a, p3b] = boundedline(NRSPHB_data.t_conv(1).t,zeros(size(NRSPHB_data.t_conv(1).t)),NRSPHB_trends_server.EEMD_conf_std_limit(1,:))
set(p3b,'FaceColor','r','FaceAlpha',0.2)
set(p3a,'LineStyle','None')

set(gca,'LineWidth',2,'Box','On','FontSize',18,'YLim',[-3 4],'XLim',[datenum(1950,01,01) datenum(2021,01,01)],...
    'XTick',[datenum(1950,01,01) datenum(1960,01,01) datenum(1970,01,01) datenum(1980,01,01) datenum(1990,01,01) ...
    datenum(2000,01,01) datenum(2010,01,01) datenum(2020,01,01)],'XTickLabels',[{'1950'} {'1960'} {'1970'} {'1980'} {'1990'} {'2000'} {'2010'} {'2020'}])

ylabel('Temperature Anomaly [^\circC]');

leg = legend([p1 p2 p3 p3b],'De-seasoned temperatures','Insignificant EEMD Trend','Significant EEMD Trend','Boundary of insignificance');
set(leg,'Location','NorthWest','Box','Off','Position',[0.101466049382716 0.767275373893183 0.238859958518986 0.146990744076638]);

axes('Parent',gcf,...
    'Position',[0.620804398148148 0.11 0.284195601851853 0.815]);

p1 = plot(NRSPHB_trends.EEMD_t_conv(1).t,NRSPHB_trends.EEMD_T{1}-tr(1),'LineWidth',2)
hold on;

imf = NRSPHB_trends.EEMD_imfs.IMF_1;
a = 0;
for n = 1:size(imf,1)
    if n < 2
        a = a+4;
        p2 = plot(NRSPHB_trends.EEMD_t_conv(1).t,imf(n,:)-a,'LineWidth',2,'Color','k')
    else
        a = a+2;
        if n < size(imf,1)
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

leg = legend([p1 p2 p3],'De-seasoned temperatures','IMFs','Trend');
set(leg,'Location','SouthWest','Box','Off');

print(gcf, '-dpng','-r400', [options.plot_dir,'EEMD_example_NRSPHB_2m'])

%% Example ITA plot

%% create the timeseries 

% % random noise
% rand_noise = cumsum(-1 + 2*round(rand(1000,1)),1);
% rand_noise = (rand_noise-min(rand_noise))/(max(rand_noise)-min(rand_noise));
% rand_noise = rand_noise*0.2;
% % Linear trend
% linear_trend = (0.05*1:1:1000)/1000;
% % exponential trend
% exp_trend = exp(1:1:5)/nanmax(exp(1:1:5))*2;
% exp_trend = interp1(1:1:5,exp_trend,1:0.0039999999999999:5,'Linear');
% exp_trend = exp_trend(1:end-1);
% % sin wave
% xsin = sin(1:10);
% xsin = interp1(1:10,xsin,1:0.009:10);
% xsin = xsin(1:end-1);
% % trends for ITA
% lin_1 = (rand_noise'+linear_trend) - nanmean(rand_noise'+linear_trend);
% lin_2 = lin_1*-1;
% e_1 = (rand_noise'+exp_trend)-nanmean(rand_noise'+exp_trend) +1;
% e_2 = e_1*-1;
% experiment = xsin+rand_noise';
% 
% figure('units','normalized','position',[0 0.1 .9 .7]);
% 
% cmap = cmocean('dense',5)
% 
% % Create axes
% axes('Parent',gcf,...
%     'Position',[0.13 0.11 0.382008101851852 0.814999999999993]);
% 
% p1 = plot(lin_1,'LineWidth',2,'Color',cmap(2,:));
% hold on;
% p1 = plot(lin_2,'LineWidth',2,'Color',cmap(3,:));
% p1 = plot(e_1,'LineWidth',2,'Color',cmap(4,:));
% p1 = plot(e_2,'LineWidth',2,'Color',cmap(5,:));
% 
% % Create axes
% axes('Parent',gcf,...
%     'Position',[0.528790509259259 0.11 0.376209490740741 0.81499999999999]);
% 
% hold on;
% % lin_1
% 
% s1 = scatter(sort(lin_1(1:500)),sort(lin_1(501:end)),'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor',cmap(2,:))
% s2 = scatter(sort(lin_2(1:500)),sort(lin_2(501:end)),'MarkerFaceColor',cmap(3,:),'MarkerEdgeColor',cmap(3,:))
% s3 = scatter(sort(e_1(1:500)),sort(e_1(501:end)),'MarkerFaceColor',cmap(4,:),'MarkerEdgeColor',cmap(4,:))
% s4 = scatter(sort(e_2(1:500)),sort(e_2(501:end)),'MarkerFaceColor',cmap(5,:),'MarkerEdgeColor',cmap(5,:))
% s5 = plot(-4:0.1:4,-4:0.1:4,'LineWidth',2,'Color','k');

%% trend plots at CH100 and BMP120

%% Heatmaps for CH100 and BMP120 ITA and MK results

% get EEMD trends

multiplier = 3653;
for n_depth = 1:11
    tr = CH100_trends.EEMD_trend{n_depth}; 
    tt = datenum(cell2mat(CH100_trends.EEMD_t(n_depth)));
    tr_rate = diff(tr);
    % 2010s
    trend_ave_CH100(n_depth).t2010s = nanmean(tr_rate)*multiplier;        
end
for n_depth = 1:12
    tr = BMP120_trends.EEMD_trend{n_depth}; 
    tt = datenum(cell2mat(BMP120_trends.EEMD_t(n_depth)));
    tr_rate = diff(tr);
    % 2010s
    trend_ave_BMP120(n_depth).t2010s = nanmean(tr_rate)*multiplier;        
end


figure('units','normalized','position',[0 0.1 .6 .7]);

cmap = cmocean('balance',23);
cmap_limits = linspace(-1.1,1.1,23);

% CH100
axes(gcf,'Position',[0.132595486111111 0.233796296296296 0.33206360479798 0.691203703703695])

ylim([1 12])
xlim([0 3])

% EEMD
for n_depth = 1:11

    tr = trend_ave_CH100(n_depth).t2010s;
    c(1) = interp1(cmap_limits,cmap(:,1),tr);
    c(2) = interp1(cmap_limits,cmap(:,2),tr);
    c(3) = interp1(cmap_limits,cmap(:,3),tr);
    color = [c(1) c(2) c(3)];
    patch([0 0 1 1],[n_depth n_depth+1 n_depth+1 n_depth],color);
    
    if n_depth == 1
        bott = 0.85;
    else
        bott = bott-0.063;
    end
    
    annotation(gcf,'textbox',...
    [0.16 bott 0.123131944444445 0.0714285714285705],...
    'String',[num2str(round(tr,2)),' *'],...
    'LineStyle','none',...
    'FontSize',14,...
    'FitBoxToText','off');

end


% ITA 
for n_depth = 1:11

    tr = CH100_trends.ITA_trend_per_decade(n_depth);
    c(1) = interp1(cmap_limits,cmap(:,1),tr);
    c(2) = interp1(cmap_limits,cmap(:,2),tr);
    c(3) = interp1(cmap_limits,cmap(:,3),tr);
    color = [c(1) c(2) c(3)];
    patch([1 1 2 2 ],[n_depth n_depth+1 n_depth+1 n_depth],color);
    
    if n_depth == 1
        bott = 0.85;
    else
        bott = bott-0.063;
    end
    
    annotation(gcf,'textbox',...
    [0.27 bott 0.123131944444445 0.0714285714285705],...
    'String',[num2str(round(tr,2)),' *'],...
    'LineStyle','none',...
    'FontSize',14,...
    'FitBoxToText','off');

end

% MK test
for n_depth = 1:11

    tr = CH100_trends.MK_trend_per_decade(n_depth);
    pv = CH100_trends.MK_pval(n_depth);
    c(1) = interp1(cmap_limits,cmap(:,1),tr);
    c(2) = interp1(cmap_limits,cmap(:,2),tr);
    c(3) = interp1(cmap_limits,cmap(:,3),tr);
    color = [c(1) c(2) c(3)];
    patch([2 2 3 3],[n_depth n_depth+1 n_depth+1 n_depth],color);
    
    if n_depth == 1
        bott = 0.85;
    else
        bott = bott-0.063;
    end
    
    if pv < 0.05
        annotation(gcf,'textbox',...
        [0.38 bott 0.123131944444445 0.0714285714285705],...
        'String',[num2str(round(tr,2)),' *'],...
        'LineStyle','none',...
        'FontSize',14,...
        'FitBoxToText','off');
    else
        annotation(gcf,'textbox',...
        [0.38 bott 0.123131944444445 0.0714285714285705],...
        'String',num2str(round(tr,2)),...
        'LineStyle','none',...
        'FontSize',14,...
        'FitBoxToText','off');
    end

end

set(gca,'YDir','Reverse','LineWidth',1,'YTick',1.5:1:11.5,'YTickLabels',[10.5, 20, 27.5, 35.5, 43.5, 51.5, 59.5, 67.5, 75.5, 84.5, 91.5],...
    'XTick',[0.5 1.5 2.5],'XTickLabels',[{'EEMD'} {'ITA'} {'TSSE'}],'FontSize',16,'Box','On');
ylabel('Depth [m]');
title('CH100')


% BMP120
axes(gcf,'Position',[0.570340909090909 0.11 0.334659090909091 0.773267195767196]);

ylim([1 13])
xlim([0 3])

% EEMD
for n_depth = 1:12

    tr = trend_ave_BMP120(n_depth).t2010s;
    c(1) = interp1(cmap_limits,cmap(:,1),tr);
    c(2) = interp1(cmap_limits,cmap(:,2),tr);
    c(3) = interp1(cmap_limits,cmap(:,3),tr);
    color = [c(1) c(2) c(3)];
    patch([0 0 1 1],[n_depth n_depth+1 n_depth+1 n_depth],color);
    
    if n_depth == 1
        bott = 0.8;
    else
        bott = bott-0.063;
    end
    
    annotation(gcf,'textbox',...
    [0.6 bott 0.123131944444445 0.0714285714285705],...
    'String',[num2str(round(tr,2)),' *'],...
    'LineStyle','none',...
    'FontSize',14,...
    'FitBoxToText','off');

end


% ITA 
for n_depth = 1:12

    tr = BMP120_trends.ITA_trend_per_decade(n_depth);
    c(1) = interp1(cmap_limits,cmap(:,1),tr);
    c(2) = interp1(cmap_limits,cmap(:,2),tr);
    c(3) = interp1(cmap_limits,cmap(:,3),tr);
    color = [c(1) c(2) c(3)];
    patch([1 1 2 2],[n_depth n_depth+1 n_depth+1 n_depth],color);
    
    if n_depth == 1
        bott = 0.8;
    else
        bott = bott-0.064;
    end
    
    annotation(gcf,'textbox',...
    [0.71 bott 0.123131944444445 0.0714285714285705],...
    'String',[num2str(round(tr,2)),' *'],...
    'LineStyle','none',...
    'FontSize',14,...
    'FitBoxToText','off');

end

% MK test
for n_depth = 1:12

    tr = BMP120_trends.MK_trend_per_decade(n_depth);
    pv = BMP120_trends.MK_pval(n_depth);
    c(1) = interp1(cmap_limits,cmap(:,1),tr);
    c(2) = interp1(cmap_limits,cmap(:,2),tr);
    c(3) = interp1(cmap_limits,cmap(:,3),tr);
    color = [c(1) c(2) c(3)];
    patch([2 2 3 3],[n_depth n_depth+1 n_depth+1 n_depth],color);
    
    if n_depth == 1
        bott = 0.8;
    else
        bott = bott-0.064;
    end
    
    if pv < 0.05
        annotation(gcf,'textbox',...
        [0.82 bott 0.123131944444445 0.0714285714285705],...
        'String',[num2str(round(tr,2)),' *'],...
        'LineStyle','none',...
        'FontSize',14,...
        'FitBoxToText','off');
    else
        annotation(gcf,'textbox',...
        [0.82 bott 0.123131944444445 0.0714285714285705],...
        'String',num2str(round(tr,2)),...
        'LineStyle','none',...
        'FontSize',14,...
        'FitBoxToText','off');
    end

end

set(gca,'YDir','Reverse','LineWidth',1,'YTick',1.5:1:12.5,'YTickLabels',[18.5, 27.5, 35, 43, 50.5, 58.5, 67, 75, 83.5, 91.5, 99.5, 107.5],...
    'XTick',[0.5 1.5 2.5],'XTickLabels',[{'EEMD'} {'ITA'} {'TSSE'}],'FontSize',16,'Box','On');
title('BMP120')

print(gcf, '-dpng','-r400', [options.plot_dir,'Trend_comparison_CH100_BMP120'])

%% Plot ITA for CH100 and BMP120

figure('units','normalized','position',[0 0.1 .8 .7]);

% Create axes
axes(gcf,'Position',[0.0818684895833333 0.11 0.397786458333333 0.814999999999983]);

cmap = cbrewer('qual','Paired',12);

clear s
for n = 1:11
  s(n) = scatter(CH100_trends.ITA_stats{n}.TEMP_half_1, CH100_trends.ITA_stats{n}.TEMP_half_2,10,'MarkerEdgeColor',cmap(n,:),'MarkerFaceColor',cmap(n,:))
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
  s(n) = scatter(BMP120_trends.ITA_stats{n}.TEMP_half_1, BMP120_trends.ITA_stats{n}.TEMP_half_2,10,'MarkerEdgeColor',cmap(n,:),'MarkerFaceColor',cmap(n,:))
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

print(gcf, '-dpng','-r400', [options.plot_dir,'ITA_comparison'])

%% Port Hacking heatmap

% get trend per decade
% MAKE SURE TO LEAVE 5 YEARS EITHER SIDE OF PERIOD FOR EDGE EFFECTS

multiplier = 12*10;
for n_depth = 1:9
    
    tr = NRSPHB_trends.EEMD_trend{n_depth}; 
    tt = datenum(cell2mat(NRSPHB_trends.EEMD_t(n_depth)));
    tr_rate = diff(tr);
    % 1960s
    check = tt >= datenum(1960,01,01) & tt < datenum(1970,01,01);
    trend_ave(n_depth).t1960s = nanmean(tr_rate(check(1:end-1)))*multiplier;    
    % 1970s
    check = tt >= datenum(1970,01,01) & tt < datenum(1980,01,01);
    trend_ave(n_depth).t1970s = nanmean(tr_rate(check(1:end-1)))*multiplier;   
    % 1980s
    check = tt >= datenum(1980,01,01) & tt < datenum(1990,01,01);
    trend_ave(n_depth).t1980s = nanmean(tr_rate(check(1:end-1)))*multiplier; 
    % 1990s
    check = tt >= datenum(1990,01,01) & tt < datenum(2000,01,01);
    trend_ave(n_depth).t1990s = nanmean(tr_rate(check(1:end-1)))*multiplier;      
    % 2000s
    check = tt >= datenum(2000,01,01) & tt < datenum(2010,01,01);
    trend_ave(n_depth).t2000s = nanmean(tr_rate(check(1:end-1)))*multiplier;    
    % 2010s
    check = tt >= datenum(2010,01,01) & tt < datenum(2020,01,01);
    trend_ave(n_depth).t2010s = nanmean(tr_rate(check(1:end-1)))*multiplier;  
    % total change
    trend_ave(n_depth).total = abs(tr(end)-tr(1));
end

figure('units','normalized','position',[0 0.1 .7 .8]);

cmap = cmocean('balance',23);
cmap_limits = -0.3:0.027:0.3;
xaxis_limits = linspace(0.1,0.9,6);
yaxis_limits = linspace(0.11,0.93,6);
xlim_limits = linspace(0,7,6);
ylim_limits = linspace(1,10,6);


ylim([1 10])
xlim([0 7])

for n_depth = 1:9

    % organise trends over time
    trs = [trend_ave(n_depth).t1960s, trend_ave(n_depth).t1970s, trend_ave(n_depth).t1980s, ...
        trend_ave(n_depth).t1990s, trend_ave(n_depth).t2000s, trend_ave(n_depth).t2010s, trend_ave(n_depth).total];
    y_pos = interp1(ylim_limits,fliplr(yaxis_limits),n_depth+1.3,'Linear','extrap');
    for n = 1:numel(trs)
        if n < 7
            c(1) = interp1(cmap_limits,cmap(:,1),trs(n));
            c(2) = interp1(cmap_limits,cmap(:,2),trs(n));
            c(3) = interp1(cmap_limits,cmap(:,3),trs(n));
            color = [c(1) c(2) c(3)];
        else
            color = [0.9 0.9 0.9];            
        end
        if n == 1
            patch([0 0 1 1 ],[n_depth n_depth+1 n_depth+1 n_depth],color);
        else
            patch([n-1 n-1 n n],[n_depth n_depth+1 n_depth+1 n_depth],color);
        end
        x_pos = interp1(xlim_limits,xaxis_limits,n-0.6);
        
        annotation(gcf,'textbox',...
        [x_pos y_pos 0.1 0.1],...
        'String',[num2str(round(trs(n),2)),' *'],...
        'LineStyle','none',...
        'FontSize',14,...
        'FitBoxToText','off');
    end

end

set(gca,'YDir','Reverse','LineWidth',1,'YTick',1.5:1:9.5,'YTickLabels',[2, 19, 31, 40, 50, 59, 75, 81, 99],...
    'XTick',0.5:1:6.5,'XTickLabels',[{'1960s'} {'1970s'} {'1980s'} {'1990s'} {'2000s'} {'2010s'} {'Total'}],'FontSize',16,'Box','On');
ylabel('Depth [m]');
title('NRS Port Hacking')

print(gcf, '-dpng','-r400', [options.plot_dir,'EEMD_trends_PortHacking'])


%% Maria Island heatmap

% get trend per decade
% MAKE SURE TO LEAVE 5 YEARS EITHER SIDE OF PERIOD FOR EDGE EFFECTS

multiplier = 12*10;
clear trend_ave
for n_depth = 1:2
    
    tr = NRSMAI_trends.EEMD_trend{n_depth}; 
    tt = datenum(cell2mat(NRSMAI_trends.EEMD_t(n_depth)));
    tr_rate = diff(tr);
   % 1950s
    check = tt >= datenum(1950,01,01) & tt < datenum(1960,01,01);
    trend_ave(n_depth).t1950s = nanmean(tr_rate(check(1:end-1)))*multiplier;        
    % 1960s
    check = tt >= datenum(1960,01,01) & tt < datenum(1970,01,01);
    trend_ave(n_depth).t1960s = nanmean(tr_rate(check(1:end-1)))*multiplier;    
    % 1970s
    check = tt >= datenum(1970,01,01) & tt < datenum(1980,01,01);
    trend_ave(n_depth).t1970s = nanmean(tr_rate(check(1:end-1)))*multiplier;   
    % 1980s
    check = tt >= datenum(1980,01,01) & tt < datenum(1990,01,01);
    trend_ave(n_depth).t1980s = nanmean(tr_rate(check(1:end-1)))*multiplier; 
    % 1990s
    check = tt >= datenum(1990,01,01) & tt < datenum(2000,01,01);
    trend_ave(n_depth).t1990s = nanmean(tr_rate(check(1:end-1)))*multiplier;      
    % 2000s
    check = tt >= datenum(2000,01,01) & tt < datenum(2010,01,01);
    trend_ave(n_depth).t2000s = nanmean(tr_rate(check(1:end-1)))*multiplier;    
    % 2010s
    check = tt >= datenum(2010,01,01) & tt < datenum(2020,01,01);
    trend_ave(n_depth).t2010s = nanmean(tr_rate(check(1:end-1)))*multiplier;     
    % total change
    trend_ave(n_depth).total = abs(tr(end)-tr(1));    
end

figure('units','normalized','position',[0 0.1 .7 .4]);

cmap = cmocean('balance',23);
cmap_limits = -0.3:0.027:0.3;
xaxis_limits = linspace(0.1,0.9,7);
yaxis_limits = linspace(0.11,0.93,2);
xlim_limits = linspace(0,8,7);
ylim_limits = linspace(1,3,2);

axes1 = axes('Parent',gcf,...
    'Position',[0.13 0.256365740740741 0.775 0.49537037037037]);


ylim([1 3])
xlim([0 8])

% ITA 
for n_depth = 1:2

    % organise trends over time
    trs = [trend_ave(n_depth).t1950s, trend_ave(n_depth).t1960s, trend_ave(n_depth).t1970s, trend_ave(n_depth).t1980s, ...
        trend_ave(n_depth).t1990s, trend_ave(n_depth).t2000s, trend_ave(n_depth).t2010s, trend_ave(n_depth).total];
    if n_depth == 1
        y_pos = 0.581;
    else
        y_pos = 0.322;
    end    
    for n = 1:numel(trs)
      if n < 8
            c(1) = interp1(cmap_limits,cmap(:,1),trs(n));
            c(2) = interp1(cmap_limits,cmap(:,2),trs(n));
            c(3) = interp1(cmap_limits,cmap(:,3),trs(n));
            color = [c(1) c(2) c(3)];
        else
            color = [0.9 0.9 0.9];            
        end
        if n == 1
            patch([0 0 1 1 ],[n_depth n_depth+1 n_depth+1 n_depth],color);
        else
            patch([n-1 n-1 n n],[n_depth n_depth+1 n_depth+1 n_depth],color);
        end
        x_pos = interp1(xlim_limits,xaxis_limits,n-0.6);
        
        annotation(gcf,'textbox',...
        [x_pos y_pos 0.1 0.1],...
        'String',[num2str(round(trs(n),2)),' *'],...
        'LineStyle','none',...
        'FontSize',14,...
        'FitBoxToText','off');
    end

end

set(gca,'YDir','Reverse','LineWidth',1,'YTick',1.5:2.5,'YTickLabels',[2, 21],...
    'XTick',0.5:1:7.5,'XTickLabels',[{'1950s'} {'1960s'} {'1970s'} {'1980s'} {'1990s'} {'2000s'} {'2010s'} {'Total'}],'FontSize',16,'Box','On');
ylabel('Depth [m]');
title('NRS Maria Island')

print(gcf, '-dpng','-r400', [options.plot_dir,'EEMD_trends_Maria Island'])

