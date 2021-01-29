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

% NRSMAI
% time
nt = size(NRSMAI_data.t,2);
for nn = 1:7
    for t = 1:nt
        a = squeeze(vertcat(NRSMAI_data.t(nn,t,:)))';
        NRSMAI_data.t_conv(nn).t(t) = datenum(convertCharsToStrings(a));
    end
end
% EEMD time
for nn = 1:7
    a = NRSMAI_trends.EEMD_t{nn};
    for t = 1:size(a,1)
        b = a(t,:);
        NRSMAI_trends.EEMD_t_conv(nn).t(t) = datenum(convertCharsToStrings(b));
    end
end

%% Fixing depths where IMFs not capturing trend properly

a = NRSMAI_trends.EEMD_imfs.IMF_1;
NRSMAI_trends.EEMD_trend{1} = a(end,:) + a(end-1,:);
NRSMAI_trends.EEMD_trend_EAC{1} = a(end,:) + a(end-1,:) + a(end-2,:);

%% calculate average trends

multiplier = 12*10; % 10 years in months

[trends_PHB] = get_trends(NRSPHB_data, NRSPHB_data_server, NRSPHB_trends,NRSPHB_trends_server,multiplier);
[trends_MAI] = get_trends(NRSMAI_data, NRSMAI_data_server, NRSMAI_trends,NRSMAI_trends_server,multiplier);

%% Figure Port Hacking

figure('units','normalized','position',[0 0.1 .7 .8]);

cmap = cmocean('balance',23);
cmap_limits = linspace(-1,1,23);
xaxis_limits = linspace(0.1,0.9,6);
yaxis_limits = linspace(0.11,0.93,6);
xlim_limits = linspace(0,6,6);
ylim_limits = linspace(1,10,6);


ylim([1 10])
xlim([0 6])

for n_depth = 1:9

    % organise trends over time
    trs_sig = [trends_PHB(n_depth).t1960s_sig, trends_PHB(n_depth).t1970s_sig, trends_PHB(n_depth).t1980s_sig, ...
        trends_PHB(n_depth).t1990s_sig, trends_PHB(n_depth).t2000s_sig, trends_PHB(n_depth).t2010s_sig, trends_PHB(n_depth).total_sig_1990s];
    trs_insig = [trends_PHB(n_depth).t1960s_insig, trends_PHB(n_depth).t1970s_insig, trends_PHB(n_depth).t1980s_insig, ...
        trends_PHB(n_depth).t1990s_insig, trends_PHB(n_depth).t2000s_insig, trends_PHB(n_depth).t2010s_insig, trends_PHB(n_depth).total_insig];
    trs_EAC_sig = [trends_PHB(n_depth).t1960s_EAC_sig, trends_PHB(n_depth).t1970s_EAC_sig, trends_PHB(n_depth).t1980s_EAC_sig, ...
        trends_PHB(n_depth).t1990s_EAC_sig, trends_PHB(n_depth).t2000s_EAC_sig, trends_PHB(n_depth).t2010s_EAC_sig, trends_PHB(n_depth).total_EAC_sig_1990s];    
    trs_EAC_insig = [trends_PHB(n_depth).t1960s_EAC_insig, trends_PHB(n_depth).t1970s_EAC_insig, trends_PHB(n_depth).t1980s_EAC_insig, ...
        trends_PHB(n_depth).t1990s_EAC_insig, trends_PHB(n_depth).t2000s_EAC_insig, trends_PHB(n_depth).t2010s_EAC_insig, trends_PHB(n_depth).total_EAC_insig];    
    y_pos = interp1(ylim_limits,fliplr(yaxis_limits),n_depth+1.1,'Linear','extrap');
    y_pos_EAC = interp1(ylim_limits,fliplr(yaxis_limits),n_depth+1.4,'Linear','extrap');
    for n = 1:numel(trs_sig)-1
        if n < 7
            % determine value used
            if isfinite(trs_sig(n))
                val = trs_sig(n);
                sig_indicator = 1;
            else
                val = trs_insig(n);
                sig_indicator = 0;
            end
            if isfinite(trs_EAC_sig(n))
                val_EAC = trs_EAC_sig(n);
                sig_indicator_EAC = 1;
            else
                val_EAC = trs_EAC_insig(n);
                sig_indicator_EAC = 0;
            end            
            c(1) = interp1(cmap_limits,cmap(:,1),val_EAC);
            c(2) = interp1(cmap_limits,cmap(:,2),val_EAC);
            c(3) = interp1(cmap_limits,cmap(:,3),val_EAC);
            color = [c(1) c(2) c(3)];
        else
            val = trs_sig(7);
            val_EAC = trs_EAC_sig(7);
            color = [1 1 1];            
        end
        if n == 1
            patch([0 0 1 1 ],[n_depth n_depth+1 n_depth+1 n_depth],color);
        else
            patch([n-1 n-1 n n],[n_depth n_depth+1 n_depth+1 n_depth],color);
        end
        x_pos = interp1(xlim_limits,xaxis_limits,n-0.6);
        
        if isnan(val)
            val = 0;
        end
        if isnan(val_EAC)
            val_EAC = 0;
        end
        
        if sig_indicator == 1
            annotation(gcf,'textbox',...
            [x_pos y_pos_EAC 0.1 0.1],...
            'String',[num2str(round(val,2)),' \color{black}{\diamondsuit}'],...
            'LineStyle','none',...
            'FontSize',14,...
            'FitBoxToText','off');
        else
            annotation(gcf,'textbox',...
            [x_pos y_pos_EAC 0.1 0.1],...
            'String',[num2str(round(val,2))],...
            'LineStyle','none',...
            'FontSize',14,...
            'FitBoxToText','off');                  
        end
        if sig_indicator_EAC == 1        
            annotation(gcf,'textbox',...
            [x_pos y_pos 0.1 0.1],...
            'String',['\bf', num2str(round(val_EAC,2)),' \color{black}{\diamondsuit}'],...
            'LineStyle','none',...
            'FontSize',14,...
            'FitBoxToText','off','FontName','Helvetica-Narrow');        
        else
            annotation(gcf,'textbox',...
            [x_pos y_pos 0.1 0.1],...
            'String',['\bf', num2str(round(val_EAC,2))],...
            'LineStyle','none',...
            'FontSize',14,...
            'FitBoxToText','off','FontName','Helvetica-Narrow');     
        end
    end

end

set(gca,'YDir','Reverse','LineWidth',1,'YTick',1.5:1:9.5,'YTickLabels',[2, 19, 31, 40, 50, 59, 75, 81, 99],...
    'XTick',0.5:1:6.5,'XTickLabels',[{'1960s'} {'1970s'} {'1980s'} {'1990s'} {'2000s'} {'2010s'} {'1990-2020'}],'FontSize',16,'Box','On');
ylabel('Depth [m]');
title('NRS Port Hacking')

print(gcf, '-dpng','-r400', [options.plot_dir,'EEMD_trends_PortHacking'])


%% Figure Maria Island

figure('units','normalized','position',[0 0.1 .7 .8]);

cmap = cmocean('balance',23);
cmap_limits = linspace(-1,1,23);
xaxis_limits = linspace(0.06,0.9,6);
yaxis_limits = linspace(0.11,0.93,6);
xlim_limits = linspace(0,8,6);
ylim_limits = linspace(1,8,6);


ylim([1 7])
xlim([0 7])
y_spacing = linspace(1.1,2,6);
x_spacing = linspace(0,0.2,6);
for n_depth = 1:6

    % organise trends over time
    trs_sig = [trends_MAI(n_depth).t1950s_sig, trends_MAI(n_depth).t1960s_sig, trends_MAI(n_depth).t1970s_sig, trends_MAI(n_depth).t1980s_sig, ...
        trends_MAI(n_depth).t1990s_sig, trends_MAI(n_depth).t2000s_sig, trends_MAI(n_depth).t2010s_sig, trends_MAI(n_depth).total_sig_1990s];
    trs_insig = [trends_MAI(n_depth).t1950s_insig, trends_MAI(n_depth).t1960s_insig, trends_MAI(n_depth).t1970s_insig, trends_MAI(n_depth).t1980s_insig, ...
        trends_MAI(n_depth).t1990s_insig, trends_MAI(n_depth).t2000s_insig, trends_MAI(n_depth).t2010s_insig, trends_MAI(n_depth).total_insig];
    trs_EAC_sig = [trends_MAI(n_depth).t1950s_EAC_sig, trends_MAI(n_depth).t1960s_EAC_sig, trends_MAI(n_depth).t1970s_EAC_sig, trends_MAI(n_depth).t1980s_EAC_sig, ...
        trends_MAI(n_depth).t1990s_EAC_sig, trends_MAI(n_depth).t2000s_EAC_sig, trends_MAI(n_depth).t2010s_EAC_sig, trends_MAI(n_depth).total_EAC_sig_1990s];    
    trs_EAC_insig = [trends_MAI(n_depth).t1950s_EAC_insig, trends_MAI(n_depth).t1960s_EAC_insig, trends_MAI(n_depth).t1970s_EAC_insig, trends_MAI(n_depth).t1980s_EAC_insig, ...
        trends_MAI(n_depth).t1990s_EAC_insig, trends_MAI(n_depth).t2000s_EAC_insig, trends_MAI(n_depth).t2010s_EAC_insig, trends_MAI(n_depth).total_EAC_insig];    
    y_pos = interp1(ylim_limits,fliplr(yaxis_limits),n_depth+y_spacing(n_depth),'Linear','extrap');
    y_pos_EAC = interp1(ylim_limits,fliplr(yaxis_limits),n_depth+y_spacing(n_depth)+0.3,'Linear','extrap');
    xx = x_spacing(n_depth);
    for n = 1:numel(trs_sig)-1
        if n <= 7
            % determine value used
            if isfinite(trs_sig(n))
                val = trs_sig(n);
                sig_indicator = 1;
            else
                val = trs_insig(n);
                sig_indicator = 0;
            end
            if isfinite(trs_EAC_sig(n))
                val_EAC = trs_EAC_sig(n);
                sig_indicator_EAC = 1;
            else
                val_EAC = trs_EAC_insig(n);
                sig_indicator_EAC = 0;
            end            
            c(1) = interp1(cmap_limits,cmap(:,1),val_EAC);
            c(2) = interp1(cmap_limits,cmap(:,2),val_EAC);
            c(3) = interp1(cmap_limits,cmap(:,3),val_EAC);
            color = [c(1) c(2) c(3)];
            if sum(isfinite(color)) == 0
                color = [1 1 1];
            end
        else
            val = trs_sig(7);
            val_EAC = trs_EAC_sig(7);
            color = [1 1 1];            
        end
        if n == 1
            patch([0 0 1 1 ],[n_depth n_depth+1 n_depth+1 n_depth],color);
        else
            patch([n-1 n-1 n n],[n_depth n_depth+1 n_depth+1 n_depth],color);
        end
        x_pos = interp1(xlim_limits,xaxis_limits,n+xx);
        
        if isnan(val)
            val = 0;
        end
        if isnan(val_EAC)
            val_EAC = 0;
        end
        
        if sig_indicator == 1
            annotation(gcf,'textbox',...
            [x_pos y_pos_EAC 0.1 0.1],...
            'String',[num2str(round(val,2)),' \color{black}{\diamondsuit}'],...
            'LineStyle','none',...
            'FontSize',14,...
            'FitBoxToText','off');
        else
            annotation(gcf,'textbox',...
            [x_pos y_pos_EAC 0.1 0.1],...
            'String',[num2str(round(val,2))],...
            'LineStyle','none',...
            'FontSize',14,...
            'FitBoxToText','off');                  
        end
        if sig_indicator_EAC == 1        
            annotation(gcf,'textbox',...
            [x_pos y_pos 0.1 0.1],...
            'String',['\bf', num2str(round(val_EAC,2)),' \color{black}{\diamondsuit}'],...
            'LineStyle','none',...
            'FontSize',14,...
            'FitBoxToText','off','FontName','Helvetica-Narrow');        
        else
            annotation(gcf,'textbox',...
            [x_pos y_pos 0.1 0.1],...
            'String',['\bf', num2str(round(val_EAC,2))],...
            'LineStyle','none',...
            'FontSize',14,...
            'FitBoxToText','off','FontName','Helvetica-Narrow');     
        end
    end

end

set(gca,'YDir','Reverse','LineWidth',1,'YTick',1.5:1:9.5,'YTickLabels',[2, 10, 20, 30, 40, 50, 85],...
    'XTick',0.5:1:6.5,'XTickLabels',[{'1950s'} {'1960s'} {'1970s'} {'1980s'} {'1990s'} {'2000s'} {'2010s'}],'FontSize',16,'Box','On');
ylabel('Depth [m]');
title('NRS Maria Island')

print(gcf, '-dpng','-r400', [options.plot_dir,'EEMD_trends_MariaIsland'])





