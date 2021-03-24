%% Load in data

addpath(genpath('C:\Users\mphem\Documents\Work\UNSW\'))
options.data_dir = 'C:\Users\mphem\Documents\Work\UNSW\Trends\Data\';
options.plot_dir = 'C:\Users\mphem\Documents\Work\UNSW\Trends\Plots\';

% NRSPHB
NRSPHB_data = load([options.data_dir,'Split_data_server_deseason\NRSPHB_data_server1.mat']);
[NRSPHB_data, NRSPHB_trends] = combine_files('NRSPHB',NRSPHB_data);
% NRSMAI
NRSMAI_data = load([options.data_dir,'Split_data_server_deseason\NRSMAI_data_server1.mat']);
[NRSMAI_data, NRSMAI_trends] = combine_files('NRSMAI',NRSMAI_data);
%CH100
CH100_data = load([options.data_dir,'Split_data_server_deseason\CH100_data_server1.mat']);
[CH100_data, CH100_trends] = combine_files('CH100',CH100_data);
%BMP120
BMP120_data = load([options.data_dir,'Split_data_server_deseason\BMP120_data_server1.mat']);
[BMP120_data, BMP120_trends] = combine_files('BMP120',BMP120_data);
% NRSNSI
NRSNSI_data = load([options.data_dir,'Split_data_server_deseason\NRSNSI_data_server1.mat']);
[NRSNSI_data, NRSNSI_trends] = combine_files('NRSNSI',NRSNSI_data);

%% Manually fixing some trends when algorithm did not work correctly

%------------------------------------------------------------------------------------
% All seasons
%--------------------
% NRSMAI depth 3
a = NRSMAI_trends.EEMD_imfs.IMF_3;
NRSMAI_trends.EEMD_trend(3,:) = a(end,:)+a(end-1,:);
NRSMAI_trends.EEMD_trend_EAC(3,:) = a(end,:)+a(end-1,:)+a(end-2,:);
%--------------------
% CH100 depth 1
a = CH100_trends.EEMD_imfs.IMF_1;
CH100_trends.EEMD_trend{1} = a(end,:);
CH100_trends.EEMD_trend_EAC{1} = a(end,:)+a(end-1,:);
%--------------------
% CH100 depth 4
a = CH100_trends.EEMD_imfs.IMF_4;
CH100_trends.EEMD_trend{4} = a(end,:);
CH100_trends.EEMD_trend_EAC{4} = a(end,:)+a(end-1,:);
%--------------------
% CH100 depth 3
a = CH100_trends.EEMD_imfs.IMF_3;
CH100_trends.EEMD_trend{3} = a(end,:);
CH100_trends.EEMD_trend_EAC{3} = a(end,:)+a(end-1,:);
%--------------------
% CH100 depth 8
a = CH100_trends.EEMD_imfs.IMF_8;
CH100_trends.EEMD_trend{8} = a(end,:);
CH100_trends.EEMD_trend_EAC{8} = a(end,:)+a(end-1,:);
%--------------------
% CH100 depth 9
a = CH100_trends.EEMD_imfs.IMF_9;
CH100_trends.EEMD_trend{9} = a(end,:)+a(end-1,:);
CH100_trends.EEMD_trend_EAC{9} = a(end,:)+a(end-1,:)+a(end-2,:);
%--------------------
% CH100 depth 11
a = CH100_trends.EEMD_imfs.IMF_11;
CH100_trends.EEMD_trend{11} = a(end,:);
CH100_trends.EEMD_trend_EAC{11} = a(end,:)+a(end-1,:);
%--------------------
% BMP120 all depths
% for n = 1:13
%     a = eval(['BMP120_trends.EEMD_imfs.IMF_',num2str(n)]);
%     BMP120_trends.EEMD_trend{n} = a(end,:)+a(end-1,:);
%     BMP120_trends.EEMD_trend_EAC{n} = a(end,:)+a(end-1,:)+a(end-2,:);
% end
%--------------------
% BMP120 depth 3
a = BMP120_trends.EEMD_imfs.IMF_3;
BMP120_trends.EEMD_trend{3} = a(end,:)+a(end-1,:);
BMP120_trends.EEMD_trend_EAC{3} = a(end,:)+a(end-1,:)+a(end-2,:);



%% calculate average trends

multiplier = 12*10; % 10 years in months
[trends_PHB] = get_trends(NRSPHB_data, NRSPHB_data, NRSPHB_trends,NRSPHB_trends,multiplier);
[trends_MAI] = get_trends(NRSMAI_data, NRSMAI_data, NRSMAI_trends,NRSMAI_trends,multiplier);
multiplier = 10*365.25; % 10 years in days
[trends_CH100] = get_trends(CH100_data, CH100_data, CH100_trends, CH100_trends,multiplier);
[trends_BMP120] = get_trends(BMP120_data, BMP120_data, BMP120_trends, BMP120_trends,multiplier);
[trends_NRSNSI] = get_trends(NRSNSI_data, NRSNSI_data, NRSNSI_trends, NRSNSI_trends,multiplier);

%% Figure Port Hacking
cmap = cmocean('balance',23);
cmap_limits = linspace(-0.9,0.9,23);
xaxis_limits = linspace(0.1,0.9,6);
yaxis_limits = linspace(0.11,0.93,6);
xlim_limits = linspace(0,6,6);
ylim_limits = linspace(1,10,6);

figure('units','normalized','position',[0 -.4 0.7 1.4]);
ylim([1 12])
xlim([0 7])

row_n = 0;
for n_depth = 1:7

    % organise trends over time
    % trends
    trs_sig = [trends_PHB(n_depth).t1960s_sig, trends_PHB(n_depth).t1970s_sig, trends_PHB(n_depth).t1980s_sig, ...
        trends_PHB(n_depth).t1990s_sig, trends_PHB(n_depth).t2000s_sig, trends_PHB(n_depth).t2010s_sig, trends_PHB(n_depth).total_sig_1990s];
    trs_insig = [trends_PHB(n_depth).t1960s_insig, trends_PHB(n_depth).t1970s_insig, trends_PHB(n_depth).t1980s_insig, ...
        trends_PHB(n_depth).t1990s_insig, trends_PHB(n_depth).t2000s_insig, trends_PHB(n_depth).t2010s_insig, trends_PHB(n_depth).total_insig];
    trs_EAC_sig = [trends_PHB(n_depth).t1960s_EAC_sig, trends_PHB(n_depth).t1970s_EAC_sig, trends_PHB(n_depth).t1980s_EAC_sig, ...
        trends_PHB(n_depth).t1990s_EAC_sig, trends_PHB(n_depth).t2000s_EAC_sig, trends_PHB(n_depth).t2010s_EAC_sig, trends_PHB(n_depth).total_EAC_sig_1990s];    
    trs_EAC_insig = [trends_PHB(n_depth).t1960s_EAC_insig, trends_PHB(n_depth).t1970s_EAC_insig, trends_PHB(n_depth).t1980s_EAC_insig, ...
        trends_PHB(n_depth).t1990s_EAC_insig, trends_PHB(n_depth).t2000s_EAC_insig, trends_PHB(n_depth).t2010s_EAC_insig, trends_PHB(n_depth).total_EAC_insig];    
    % stds
    trs_sig_std = [trends_PHB(n_depth).t1960s_sig_std, trends_PHB(n_depth).t1970s_sig_std, trends_PHB(n_depth).t1980s_sig_std, ...
        trends_PHB(n_depth).t1990s_sig_std, trends_PHB(n_depth).t2000s_sig_std, trends_PHB(n_depth).t2010s_sig_std];
    trs_insig_std = [trends_PHB(n_depth).t1960s_insig_std, trends_PHB(n_depth).t1970s_insig_std, trends_PHB(n_depth).t1980s_insig_std, ...
        trends_PHB(n_depth).t1990s_insig_std, trends_PHB(n_depth).t2000s_insig_std, trends_PHB(n_depth).t2010s_insig_std];
    trs_EAC_sig_std = [trends_PHB(n_depth).t1960s_EAC_sig_std, trends_PHB(n_depth).t1970s_EAC_sig_std, trends_PHB(n_depth).t1980s_EAC_sig_std, ...
        trends_PHB(n_depth).t1990s_EAC_sig_std, trends_PHB(n_depth).t2000s_EAC_sig_std, trends_PHB(n_depth).t2010s_EAC_sig_std];    
    trs_EAC_insig_std = [trends_PHB(n_depth).t1960s_EAC_insig_std, trends_PHB(n_depth).t1970s_EAC_insig_std, trends_PHB(n_depth).t1980s_EAC_insig_std, ...
        trends_PHB(n_depth).t1990s_EAC_insig_std, trends_PHB(n_depth).t2000s_EAC_insig_std, trends_PHB(n_depth).t2010s_EAC_insig_std];    

    row_n = row_n+1;
    %---------------------------------------------------------------------------------------------------
    % EEMD trends
    for n = 1:numel(trs_sig)-1
        % determine value used
        if isfinite(trs_sig(n))
            val = trs_sig(n);
            std = trs_sig_std(n);
            sig_indicator = 1;
        else
            val = trs_insig(n);
            std = trs_insig_std(n);
            sig_indicator = 0;
        end
        if isfinite(trs_EAC_sig(n))
            val_EAC = trs_EAC_sig(n);
            std_EAC = trs_EAC_sig_std(n);
            sig_indicator_EAC = 1;
        else
            val_EAC = trs_EAC_insig(n);
            std_EAC = trs_EAC_insig_std(n);
            sig_indicator_EAC = 0;
        end                     
        
        c(1) = interp1(cmap_limits,cmap(:,1),val_EAC);
        c(2) = interp1(cmap_limits,cmap(:,2),val_EAC);
        c(3) = interp1(cmap_limits,cmap(:,3),val_EAC);
        color_EAC = [c(1) c(2) c(3)];        
        c(1) = interp1(cmap_limits,cmap(:,1),val);
        c(2) = interp1(cmap_limits,cmap(:,2),val);
        c(3) = interp1(cmap_limits,cmap(:,3),val);
        color = [c(1) c(2) c(3)];              
        
        if isnan(val)
            val = 0;
            std = ' ';
        end
        if isnan(val_EAC)
            val_EAC = 0;
            std_EAC = ' ';
        end
        val_str = num2str(round(val,6));
        val_EAC_str = num2str(round(val_EAC,6));
        std_str = num2str(round(std,6));
        std_EAC_str = num2str(round(std_EAC,6));        
        % ensure value is same size
        
        
        if isempty(strfind(val_str,'-')) 
            val_str = val_str(1:4);
        else
            val_str = val_str(1:5);
        end
        std_str = std_str(1:4);
        if isempty(strfind(val_EAC_str,'-')) 
            val_EAC_str = val_EAC_str(1:4);
        else
            val_EAC_str = val_EAC_str(1:5);
        end
        std_EAC_str = std_EAC_str(1:4);
        
        if ~isempty(strfind(std_str,'e'))
            std_str = '0.00';
        end
        if ~isempty(strfind(std_EAC_str,'e'))
            std_EAC_str = '0.00';
        end
        
        if sig_indicator == 1      
            text(n-0.7, row_n+0.7, [val_str,' \pm ',std_str,''],'BackgroundColor',color, ...
                'EdgeColor',[0 0 0],'LineWidth',1)
        else      
            text(n-0.7, row_n+0.7, [val_str,' \pm ',std_str,''],'BackgroundColor',color,'color',[0 0 0])
        end
        
        if sig_indicator_EAC == 1             
            text(n-0.7, row_n+0.3, ['\bf',val_EAC_str,' \pm ',std_EAC_str,''], ...
                'BackgroundColor',color_EAC,'EdgeColor',[0 0 0],'LineWidth',2)
        else                  
            text(n-0.7, row_n+0.3, ['\bf',val_EAC_str,' \pm ',std_EAC_str,''],'BackgroundColor',color_EAC,'color',[0 0 0])
        end        
    end
    %-------------------------------------------------------------------------------------------------------------------
    % Adding MK trends at the end for whole time period
    
    MK = NRSPHB_trends.MK_trend_per_decade(n_depth);
    ps = NRSPHB_trends.MK_pval(n_depth);
    % get MK color
    c(1) = interp1(cmap_limits,cmap(:,1),MK);
    c(2) = interp1(cmap_limits,cmap(:,2),MK);
    c(3) = interp1(cmap_limits,cmap(:,3),MK);
    color = [c(1) c(2) c(3)];      
    if ps < 0.05
        text(n+0.2, row_n+0.5, [num2str(round(MK,2))],'BackgroundColor',color, ...
            'EdgeColor',[0 0 0],'LineWidth',1);
    else
        text(n+0.2, row_n+0.5, [num2str(round(MK,2))],'BackgroundColor',color);        
    end

end
 
set(gca,'YDir','Reverse','Visible','Off')

% set(gca,'YDir','Reverse','LineWidth',1,'YTick',1.5:1:9.5,'YTickLabels',[2, 19, 31, 40, 50, 59, 75, 81, 99],...
%     'XTick',0.5:1:6.5,'XTickLabels',[{'1960s'} {'1970s'} {'1980s'} {'1990s'} {'2000s'} {'2010s'} {'1990-2020'}],'FontSize',16,'Box','On');
% ylabel('Depth [m]');
% title('NRS Port Hacking')

%%
print(gcf, '-dpng','-r400', [options.plot_dir,'EEMD_trends_PortHacking'])
close all

%% Figure Maria Island

figure('units','normalized','position',[0 -.4 0.7 1.4]);

cmap = cmocean('balance',23);
cmap_limits = linspace(-1.5,1.5,23);
xaxis_limits = linspace(0.06,0.9,6);
yaxis_limits = linspace(0.11,0.93,6);
xlim_limits = linspace(0,8,6);
ylim_limits = linspace(1,8,6);
ylim([1 12])
xlim([0 8])
y_spacing = linspace(1.1,2,6);
x_spacing = linspace(0,0.2,6);


row_n = 0;
for n_depth = 1:3

    % organise trends over time
    % trends
    trs_sig = [trends_MAI(n_depth).t1950s_sig, trends_MAI(n_depth).t1960s_sig, trends_MAI(n_depth).t1970s_sig, trends_MAI(n_depth).t1980s_sig, ...
        trends_MAI(n_depth).t1990s_sig, trends_MAI(n_depth).t2000s_sig, trends_MAI(n_depth).t2010s_sig, trends_MAI(n_depth).total_sig_1990s];
    trs_insig = [trends_MAI(n_depth).t1950s_insig, trends_MAI(n_depth).t1960s_insig, trends_MAI(n_depth).t1970s_insig, trends_MAI(n_depth).t1980s_insig, ...
        trends_MAI(n_depth).t1990s_insig, trends_MAI(n_depth).t2000s_insig, trends_MAI(n_depth).t2010s_insig, trends_MAI(n_depth).total_insig];
    trs_EAC_sig = [trends_MAI(n_depth).t1950s_EAC_sig, trends_MAI(n_depth).t1960s_EAC_sig, trends_MAI(n_depth).t1970s_EAC_sig, trends_MAI(n_depth).t1980s_EAC_sig, ...
        trends_MAI(n_depth).t1990s_EAC_sig, trends_MAI(n_depth).t2000s_EAC_sig, trends_MAI(n_depth).t2010s_EAC_sig, trends_MAI(n_depth).total_EAC_sig_1990s];    
    trs_EAC_insig = [trends_MAI(n_depth).t1950s_EAC_insig, trends_MAI(n_depth).t1960s_EAC_insig, trends_MAI(n_depth).t1970s_EAC_insig, trends_MAI(n_depth).t1980s_EAC_insig, ...
        trends_MAI(n_depth).t1990s_EAC_insig, trends_MAI(n_depth).t2000s_EAC_insig, trends_MAI(n_depth).t2010s_EAC_insig, trends_MAI(n_depth).total_EAC_insig];    
    % stds
    trs_sig_std = [trends_MAI(n_depth).t1950s_sig_std, trends_MAI(n_depth).t1960s_sig_std, trends_MAI(n_depth).t1970s_sig_std, trends_MAI(n_depth).t1980s_sig_std, ...
        trends_MAI(n_depth).t1990s_sig_std, trends_MAI(n_depth).t2000s_sig_std, trends_MAI(n_depth).t2010s_sig_std];
    trs_insig_std = [trends_MAI(n_depth).t1950s_insig_std, trends_MAI(n_depth).t1960s_insig_std, trends_MAI(n_depth).t1970s_insig_std, trends_MAI(n_depth).t1980s_insig_std, ...
        trends_MAI(n_depth).t1990s_insig_std, trends_MAI(n_depth).t2000s_insig_std, trends_MAI(n_depth).t2010s_insig_std];
    trs_EAC_sig_std = [trends_MAI(n_depth).t1950s_EAC_sig_std, trends_MAI(n_depth).t1960s_EAC_sig_std, trends_MAI(n_depth).t1970s_EAC_sig_std, trends_MAI(n_depth).t1980s_EAC_sig_std, ...
        trends_MAI(n_depth).t1990s_EAC_sig_std, trends_MAI(n_depth).t2000s_EAC_sig_std, trends_MAI(n_depth).t2010s_EAC_sig_std];    
    trs_EAC_insig_std = [trends_MAI(n_depth).t1950s_EAC_insig_std, trends_MAI(n_depth).t1960s_EAC_insig_std, trends_MAI(n_depth).t1970s_EAC_insig_std, trends_MAI(n_depth).t1980s_EAC_insig_std, ...
        trends_MAI(n_depth).t1990s_EAC_insig_std, trends_MAI(n_depth).t2000s_EAC_insig_std, trends_MAI(n_depth).t2010s_EAC_insig_std];    

    row_n = row_n+1;
    %-------------------------------------------------------------------------------------------------------------------
    % Adding trend  boxes    
    for n = 1:numel(trs_sig)-1
        % determine value used
        if isfinite(trs_sig(n))
            val = trs_sig(n);
            std = trs_sig_std(n);
            sig_indicator = 1;
        else
            val = trs_insig(n);
            std = trs_insig_std(n);
            sig_indicator = 0;
        end
        if isfinite(trs_EAC_sig(n))
            val_EAC = trs_EAC_sig(n);
            std_EAC = trs_EAC_sig_std(n);
            sig_indicator_EAC = 1;
        else
            val_EAC = trs_EAC_insig(n);
            std_EAC = trs_EAC_insig_std(n);
            sig_indicator_EAC = 0;
        end                     
        
        c(1) = interp1(cmap_limits,cmap(:,1),val_EAC);
        c(2) = interp1(cmap_limits,cmap(:,2),val_EAC);
        c(3) = interp1(cmap_limits,cmap(:,3),val_EAC);
        color_EAC = [c(1) c(2) c(3)];        
        c(1) = interp1(cmap_limits,cmap(:,1),val);
        c(2) = interp1(cmap_limits,cmap(:,2),val);
        c(3) = interp1(cmap_limits,cmap(:,3),val);
        color = [c(1) c(2) c(3)];              
        
        if isnan(val)
            val = 0;
            std = ' ';
        end
        if isnan(val_EAC)
            val_EAC = 0;
            std_EAC = ' ';
        end
        val_str = num2str(round(val,3));
        val_EAC_str = num2str(round(val_EAC,3));
        std_str = num2str(round(std,6));
        std_EAC_str = num2str(round(std_EAC,6));        
        % ensure value is same size
        
        
        if isempty(strfind(val_str,'-')) 
            val_str = val_str(1:4);
        else
            val_str = val_str(1:5);
        end
        std_str = std_str(1:4);
        if numel(val_EAC_str) < 4
            val_EAC_str = [val_EAC_str,'0'];
        end
        if isempty(strfind(val_EAC_str,'-')) 
            val_EAC_str = val_EAC_str(1:4);
        else
            val_EAC_str = val_EAC_str(1:5);
        end
        std_EAC_str = std_EAC_str(1:4);
        
        if ~isempty(strfind(std_str,'e'))
            std_str = '0.00';
        end
        if ~isempty(strfind(std_EAC_str,'e'))
            std_EAC_str = '0.00';
        end
        
        if sig_indicator == 1      
            text(n-0.7, row_n+0.7, [val_str,' \pm ',std_str,'  '],'BackgroundColor',color, ...
                'EdgeColor',[0 0 0],'LineWidth',1)
        else      
            text(n-0.7, row_n+0.7, [val_str,' \pm ',std_str,'  '],'BackgroundColor',color,'color',[0 0 0])
        end
        
        if sig_indicator_EAC == 1             
            text(n-0.7, row_n+0.3, ['\bf',val_EAC_str,' \pm ',std_EAC_str,'  '], ...
                'BackgroundColor',color_EAC,'EdgeColor',[0 0 0],'LineWidth',2)
        else                  
            text(n-0.7, row_n+0.3, ['\bf',val_EAC_str,' \pm ',std_EAC_str,'  '],'BackgroundColor',color_EAC,'color',[0 0 0])
        end        
    end
    %-------------------------------------------------------------------------------------------------------------------
    % Adding MK trends at the end for whole time period
    MK = NRSMAI_trends.MK_trend_per_decade(n_depth);
    ps = NRSMAI_trends.MK_pval(n_depth);
    % get MK color
    c(1) = interp1(cmap_limits,cmap(:,1),MK);
    c(2) = interp1(cmap_limits,cmap(:,2),MK);
    c(3) = interp1(cmap_limits,cmap(:,3),MK);
    color = [c(1) c(2) c(3)];      
    if ps < 0.05
        text(n+0.4, row_n+0.5, [num2str(round(MK,2))],'BackgroundColor',color, ...
            'EdgeColor',[0 0 0],'LineWidth',1);
    else
        text(n+0.4, row_n+0.5, [num2str(round(MK,2))],'BackgroundColor',color);        
    end
end
 
set(gca,'YDir','Reverse','Visible','Off');

clearvars -except *_data *_trends options trends_* n

%%
print(gcf, '-dpng','-r400', [options.plot_dir,'EEMD_trends_MariaIsland'])
close all

%% Figure CH100

figure('units','normalized','position',[0 -.4 0.7 1.4]);

cmap = cmocean('balance',23);
cmap_limits = linspace(-1.5,1.5,23);
xaxis_limits = linspace(0.06,0.9,6);
yaxis_limits = linspace(0.11,0.93,6);
xlim_limits = linspace(0,8,6);
ylim_limits = linspace(1,8,6);
ylim([1 12])
xlim([0 7])

row_n = 0;
n = 7;
for n_depth = 1:12

    % organise trends over time
    % trends
    trs_sig = [trends_CH100(n_depth).t2010s_sig];
    trs_insig = [trends_CH100(n_depth).t2010s_insig];
    trs_EAC_sig = [trends_CH100(n_depth).t2010s_EAC_sig];    
    trs_EAC_insig = [trends_CH100(n_depth).t2010s_EAC_insig];    
    % stds
    trs_sig_std = [trends_CH100(n_depth).t2010s_sig_std];
    trs_insig_std = [trends_CH100(n_depth).t2010s_insig_std];
    trs_EAC_sig_std = [trends_CH100(n_depth).t2010s_EAC_sig_std];    
    trs_EAC_insig_std = [trends_CH100(n_depth).t2010s_EAC_insig_std];    

    row_n = row_n+1;
    
    % determine value used
    if isfinite(trs_sig)
        val = trs_sig;
        std = trs_sig_std;
        sig_indicator = 1;
    else
        val = trs_insig;
        std = trs_insig_std;
        sig_indicator = 0;
    end
    if isfinite(trs_EAC_sig)
        val_EAC = trs_EAC_sig;
        std_EAC = trs_EAC_sig_std;
        sig_indicator_EAC = 1;
    else
        val_EAC = trs_EAC_insig;
        std_EAC = trs_EAC_insig_std;
        sig_indicator_EAC = 0;
    end                     

    c(1) = interp1(cmap_limits,cmap(:,1),val_EAC);
    c(2) = interp1(cmap_limits,cmap(:,2),val_EAC);
    c(3) = interp1(cmap_limits,cmap(:,3),val_EAC);
    color_EAC = [c(1) c(2) c(3)];        
    c(1) = interp1(cmap_limits,cmap(:,1),val);
    c(2) = interp1(cmap_limits,cmap(:,2),val);
    c(3) = interp1(cmap_limits,cmap(:,3),val);
    color = [c(1) c(2) c(3)];              

    if isnan(val)
        val = 0;
        std = ' ';
    end
    if isnan(val_EAC)
        val_EAC = 0;
        std_EAC = ' ';
    end
    val_str = num2str(round(val,6));
    val_EAC_str = num2str(round(val_EAC,6));
    std_str = num2str(round(std,6));
    std_EAC_str = num2str(round(std_EAC,6));        
    % ensure value is same size


    if isempty(strfind(val_str,'-')) 
        val_str = val_str(1:4);
    else
        val_str = val_str(1:5);
    end
    std_str = std_str(1:4);
    if isempty(strfind(val_EAC_str,'-')) 
        val_EAC_str = val_EAC_str(1:4);
    else
        val_EAC_str = val_EAC_str(1:5);
    end
    std_EAC_str = std_EAC_str(1:4);

    if ~isempty(strfind(std_str,'e'))
        std_str = '0.00';
    end
    if ~isempty(strfind(std_EAC_str,'e'))
        std_EAC_str = '0.00';
    end

    if sig_indicator == 1      
        text(n-0.7, row_n+0.7, [val_str,' \pm ',std_str,''],'BackgroundColor',color, ...
            'EdgeColor',[0 0 0],'LineWidth',1)
    else      
        text(1-0.7, row_n+0.7, [val_str,' \pm ',std_str,''],'BackgroundColor',color,'color',[0 0 0])
    end

    if sig_indicator_EAC == 1             
        text(n-0.7, row_n+0.3, ['\bf',val_EAC_str,' \pm ',std_EAC_str,''], ...
            'BackgroundColor',color_EAC,'EdgeColor',[0 0 0],'LineWidth',2)
    else                  
        text(1-0.7, row_n+0.3, ['\bf',val_EAC_str,' \pm ',std_EAC_str,''],'BackgroundColor',color_EAC,'color',[0 0 0])
    end        
    %-------------------------------------------------------------------------------------------------------------------
    % Adding MK trends at the end for whole time period
    MK = CH100_trends.MK_trend_per_decade(n_depth);
    ps = CH100_trends.MK_pval(n_depth);
    % get MK color
    c(1) = interp1(cmap_limits,cmap(:,1),MK);
    c(2) = interp1(cmap_limits,cmap(:,2),MK);
    c(3) = interp1(cmap_limits,cmap(:,3),MK);    
    color = [c(1) c(2) c(3)];      
    if ps < 0.05
        text(n-5.8, row_n+0.5, [num2str(round(MK,2))],'BackgroundColor',color, ...
            'EdgeColor',[0 0 0],'LineWidth',1);
    else
        text(n-5.8, row_n+0.5, [num2str(round(MK,2))],'BackgroundColor',color);        
    end

end
 
set(gca,'YDir','Reverse','Visible','Off');

clearvars -except *_data *_trends options trends_* n 

%%
print(gcf, '-dpng','-r400', [options.plot_dir,'EEMD_trends_CH100'])
close all

%% Figure BMP120

figure('units','normalized','position',[0 -.4 0.7 1.4]);

cmap = cmocean('balance',23);
cmap_limits = linspace(-1.5,1.5,23);
xaxis_limits = linspace(0.06,0.9,6);
yaxis_limits = linspace(0.11,0.93,6);
xlim_limits = linspace(0,8,6);
ylim_limits = linspace(1,8,6);
ylim([0 13])
xlim([0 7])

row_n = 0;
n = 7;
for n_depth = 1:13

    % organise trends over time
    % trends
    trs_sig = [trends_BMP120(n_depth).t2010s_sig];
    trs_insig = [trends_BMP120(n_depth).t2010s_insig];
    trs_EAC_sig = [trends_BMP120(n_depth).t2010s_EAC_sig];    
    trs_EAC_insig = [trends_BMP120(n_depth).t2010s_EAC_insig];    
    % stds
    trs_sig_std = [trends_BMP120(n_depth).t2010s_sig_std];
    trs_insig_std = [trends_BMP120(n_depth).t2010s_insig_std];
    trs_EAC_sig_std = [trends_BMP120(n_depth).t2010s_EAC_sig_std];    
    trs_EAC_insig_std = [trends_BMP120(n_depth).t2010s_EAC_insig_std];    

    row_n = row_n+1;
    
    % determine value used
    if isfinite(trs_sig)
        val = trs_sig;
        std = trs_sig_std;
        sig_indicator = 1;
    else
        val = trs_insig;
        std = trs_insig_std;
        sig_indicator = 0;
    end
    if isfinite(trs_EAC_sig)
        val_EAC = trs_EAC_sig;
        std_EAC = trs_EAC_sig_std;
        sig_indicator_EAC = 1;
    else
        val_EAC = trs_EAC_insig;
        std_EAC = trs_EAC_insig_std;
        sig_indicator_EAC = 0;
    end                     

    c(1) = interp1(cmap_limits,cmap(:,1),val_EAC);
    c(2) = interp1(cmap_limits,cmap(:,2),val_EAC);
    c(3) = interp1(cmap_limits,cmap(:,3),val_EAC);
    color_EAC = [c(1) c(2) c(3)];        
    c(1) = interp1(cmap_limits,cmap(:,1),val);
    c(2) = interp1(cmap_limits,cmap(:,2),val);
    c(3) = interp1(cmap_limits,cmap(:,3),val);
    color = [c(1) c(2) c(3)];              

    if isnan(val)
        val = 0;
        std = ' ';
    end
    if isnan(val_EAC)
        val_EAC = 0;
        std_EAC = ' ';
    end
    val_str = num2str(round(val,6));
    val_EAC_str = num2str(round(val_EAC,6));
    std_str = num2str(round(std,6));
    std_EAC_str = num2str(round(std_EAC,6));        
    % ensure value is same size


    if isempty(strfind(val_str,'-')) 
        val_str = val_str(1:4);
    else
        val_str = val_str(1:5);
    end
    std_str = std_str(1:4);
    if isempty(strfind(val_EAC_str,'-')) 
        val_EAC_str = val_EAC_str(1:4);
    else
        val_EAC_str = val_EAC_str(1:5);
    end
    std_EAC_str = std_EAC_str(1:4);

    if ~isempty(strfind(std_str,'e'))
        std_str = '0.00';
    end
    if ~isempty(strfind(std_EAC_str,'e'))
        std_EAC_str = '0.00';
    end

    if sig_indicator == 1      
        text(n-0.7, row_n+0.7, [val_str,' \pm ',std_str,''],'BackgroundColor',color, ...
            'EdgeColor',[0 0 0],'LineWidth',1)
    else      
        text(1-0.7, row_n+0.7, [val_str,' \pm ',std_str,''],'BackgroundColor',color,'color',[0 0 0])
    end

    if sig_indicator_EAC == 1             
        text(n-0.7, row_n+0.3, ['\bf',val_EAC_str,' \pm ',std_EAC_str,''], ...
            'BackgroundColor',color_EAC,'EdgeColor',[0 0 0],'LineWidth',2)
    else                  
        text(1-0.7, row_n+0.3, ['\bf',val_EAC_str,' \pm ',std_EAC_str,''],'BackgroundColor',color_EAC,'color',[0 0 0])
    end        
    %-------------------------------------------------------------------------------------------------------------------
    % Adding MK trends at the end for whole time period
    MK = BMP120_trends.MK_trend_per_decade(n_depth);
    ps = BMP120_trends.MK_pval(n_depth);
    % get MK color
    c(1) = interp1(cmap_limits,cmap(:,1),MK);
    c(2) = interp1(cmap_limits,cmap(:,2),MK);
    c(3) = interp1(cmap_limits,cmap(:,3),MK);    
    color = [c(1) c(2) c(3)];      
    if ps < 0.05
        text(n+0.2, row_n+0.5, [num2str(round(MK,2))],'BackgroundColor',color, ...
            'EdgeColor',[0 0 0],'LineWidth',1);
    else
        text(n+0.2, row_n+0.5, [num2str(round(MK,2))],'BackgroundColor',color);        
    end

end
 
set(gca,'YDir','Reverse','Visible','Off');

clearvars -except *_data *_trends options trends_* n

%%
print(gcf, '-dpng','-r400', [options.plot_dir,'EEMD_trends_BMP120'])
close all

%% Figure NRSNSI
cmap = cmocean('balance',23);
cmap_limits = linspace(-1.5,1.5,23);
xaxis_limits = linspace(0.1,0.9,6);
yaxis_limits = linspace(0.11,0.93,6);
xlim_limits = linspace(0,6,6);
ylim_limits = linspace(1,10,6);
    `   
figure('units','normalized','position',[0 -.4 0.7 1.4]);
ylim([1 12])
xlim([0 7])

row_n = 0;
n = 7;
for n_depth = 1:3

    % organise trends over time
    % trends
    trs_sig = [trends_NRSNSI(n_depth).t2010s_sig];
    trs_insig = [trends_NRSNSI(n_depth).t2010s_insig];
    trs_EAC_sig = [trends_NRSNSI(n_depth).t2010s_EAC_sig];    
    trs_EAC_insig = [trends_NRSNSI(n_depth).t2010s_EAC_insig];    
    % stds
    trs_sig_std = [trends_NRSNSI(n_depth).t2010s_sig_std];
    trs_insig_std = [trends_NRSNSI(n_depth).t2010s_insig_std];
    trs_EAC_sig_std = [trends_NRSNSI(n_depth).t2010s_EAC_sig_std];    
    trs_EAC_insig_std = [trends_NRSNSI(n_depth).t2010s_EAC_insig_std];    

    row_n = row_n+1;
    
    % determine value used
    if isfinite(trs_sig)
        val = trs_sig;
        std = trs_sig_std;
        sig_indicator = 1;
    else
        val = trs_insig;
        std = trs_insig_std;
        sig_indicator = 0;
    end
    if isfinite(trs_EAC_sig)
        val_EAC = trs_EAC_sig;
        std_EAC = trs_EAC_sig_std;
        sig_indicator_EAC = 1;
    else
        val_EAC = trs_EAC_insig;
        std_EAC = trs_EAC_insig_std;
        sig_indicator_EAC = 0;
    end                     

    c(1) = interp1(cmap_limits,cmap(:,1),val_EAC);
    c(2) = interp1(cmap_limits,cmap(:,2),val_EAC);
    c(3) = interp1(cmap_limits,cmap(:,3),val_EAC);
    color_EAC = [c(1) c(2) c(3)];        
    c(1) = interp1(cmap_limits,cmap(:,1),val);
    c(2) = interp1(cmap_limits,cmap(:,2),val);
    c(3) = interp1(cmap_limits,cmap(:,3),val);
    color = [c(1) c(2) c(3)];              

    if isnan(val)
        val = 0;
        std = ' ';
    end
    if isnan(val_EAC)
        val_EAC = 0;
        std_EAC = ' ';
    end
    val_str = num2str(round(val,6));
    val_EAC_str = num2str(round(val_EAC,6));
    std_str = num2str(round(std,6));
    std_EAC_str = num2str(round(std_EAC,6));        
    % ensure value is same size


    if isempty(strfind(val_str,'-')) 
        val_str = val_str(1:4);
    else
        val_str = val_str(1:5);
    end
    std_str = std_str(1:4);
    if isempty(strfind(val_EAC_str,'-')) 
        val_EAC_str = val_EAC_str(1:4);
    else
        val_EAC_str = val_EAC_str(1:5);
    end
    std_EAC_str = std_EAC_str(1:4);

    if ~isempty(strfind(std_str,'e'))
        std_str = '0.00';
    end
    if ~isempty(strfind(std_EAC_str,'e'))
        std_EAC_str = '0.00';
    end

    if sig_indicator == 1      
        text(n-0.7, row_n+0.7, [val_str,' \pm ',std_str,''],'BackgroundColor',color, ...
            'EdgeColor',[0 0 0],'LineWidth',1)
    else      
        text(1-0.7, row_n+0.7, [val_str,' \pm ',std_str,''],'BackgroundColor',color,'color',[0 0 0])
    end

    if sig_indicator_EAC == 1             
        text(n-0.7, row_n+0.3, ['\bf',val_EAC_str,' \pm ',std_EAC_str,''], ...
            'BackgroundColor',color_EAC,'EdgeColor',[0 0 0],'LineWidth',2)
    else                  
        text(1-0.7, row_n+0.3, ['\bf',val_EAC_str,' \pm ',std_EAC_str,''],'BackgroundColor',color_EAC,'color',[0 0 0])
    end        
    %-------------------------------------------------------------------------------------------------------------------
    % Adding MK trends at the end for whole time period
    MK = NRSNSI_trends.MK_trend_per_decade(n_depth);
    ps = NRSNSI_trends.MK_pval(n_depth);
    % get MK color
    c(1) = interp1(cmap_limits,cmap(:,1),MK);
    c(2) = interp1(cmap_limits,cmap(:,2),MK);
    c(3) = interp1(cmap_limits,cmap(:,3),MK);    
    color = [c(1) c(2) c(3)];          
    if ps < 0.05
        text(n-5.8, row_n+0.5, [num2str(round(MK,2))],'BackgroundColor',color, ...
            'EdgeColor',[0 0 0],'LineWidth',1);
    else
        text(n-5.8, row_n+0.5, [num2str(round(MK,2))],'BackgroundColor',color);        
    end

end
 
set(gca,'YDir','Reverse','Visible','Off');

clearvars -except *_data *_trends options trends_*

%%
print(gcf, '-dpng','-r400', [options.plot_dir,'EEMD_trends_NRSNSI'])
close all

%% Create Large Australia Coast


API = load('C:\Users\mphem\Documents\Work\UNSW\climatology\Revamped_scripts\Climatology\Utilities\Plot_google_map\api_key');
  
style = ['https://maps.googleapis.com/maps/api/staticmap?key=',...
    API.apiKey,'&center=-34.00179784835884,150.98879117742268&', ...
    'zoom=10&format=png&maptype=roadmap&style=element:labels%7C', ...
    'visibility:off&style=feature:administrative%7Celement:geometry%7Cvisibility:', ...
    'off&style=feature:administrative.land_parcel%7Cvisibility:off&style=feature:administrative.', ...
    'neighborhood%7Cvisibility:off&style=feature:landscape%7Ccolor:0x000000&style=feature:', ...
    'poi%7Cvisibility:off&style=feature:road%7Cvisibility:off&style=feature:road%7Celement:labels.', ...
    'icon%7Cvisibility:off&style=feature:transit%7Cvisibility:off&style=feature:water%7Ccolor:0xffffff&size=480x360'];

% the reference stations
Locations.NRSPHB = [151.216 -34.118];
Locations.NRSMAI = [148.23 -42.6];
Locations.NRSNSI = [153.562 -27.342]
Locations.BMP120 = [150.3134 -36.2064];
Locations.CH100 = [153.3957 -30.2656];


figure('units','normalized','position',[.1 0 .25 .9])

xlim([153 154])
ylim([-45 -15])
plot_google_map('Style',style,'Scale',2)
set(gca,'Visible','Off')
set(gcf,'Color','w')

hold on
scatter(Locations.NRSPHB(1),Locations.NRSPHB(2),100,'filled','MarkerFaceColor','k','MarkerEdgeColor','k','Marker','Sq')
scatter(Locations.NRSPHB(1),Locations.NRSPHB(2),50,'filled','MarkerFaceColor','w','MarkerEdgeColor','w','Marker','Sq')
scatter(Locations.NRSMAI(1),Locations.NRSMAI(2),100,'filled','MarkerFaceColor','k','MarkerEdgeColor','k','Marker','Sq')
scatter(Locations.NRSMAI(1),Locations.NRSMAI(2),50,'filled','MarkerFaceColor','w','MarkerEdgeColor','w','Marker','Sq')
scatter(Locations.BMP120(1),Locations.BMP120(2),100,'filled','MarkerFaceColor','k','MarkerEdgeColor','k','Marker','Sq')
scatter(Locations.BMP120(1),Locations.BMP120(2),50,'filled','MarkerFaceColor','w','MarkerEdgeColor','w','Marker','Sq')
scatter(Locations.CH100(1),Locations.CH100(2),100,'filled','MarkerFaceColor','k','MarkerEdgeColor','k','Marker','Sq')
scatter(Locations.CH100(1),Locations.CH100(2),50,'filled','MarkerFaceColor','w','MarkerEdgeColor','w','Marker','Sq')
scatter(Locations.NRSNSI(1),Locations.NRSNSI(2),100,'filled','MarkerFaceColor','k','MarkerEdgeColor','k','Marker','Sq')
scatter(Locations.NRSNSI(1),Locations.NRSNSI(2),50,'filled','MarkerFaceColor','w','MarkerEdgeColor','w','Marker','Sq')

%%
print(gcf, '-dpng','-r400', [options.plot_dir,'EEMD_trends_Coast'])

%% old tables

% for n_depth = 1:9
%     for n = 1:numel(trs_sig)-1
%         if n < 7
%             % determine value used
%             if isfinite(trs_sig(n))
%                 val = trs_sig(n);
%                 sig_indicator = 1;
%             else
%                 val = trs_insig(n);
%                 sig_indicator = 0;
%             end
%             if isfinite(trs_EAC_sig(n))
%                 val_EAC = trs_EAC_sig(n);
%                 sig_indicator_EAC = 1;
%             else
%                 val_EAC = trs_EAC_insig(n);
%                 sig_indicator_EAC = 0;
%             end            
%             c(1) = interp1(cmap_limits,cmap(:,1),val_EAC);
%             c(2) = interp1(cmap_limits,cmap(:,2),val_EAC);
%             c(3) = interp1(cmap_limits,cmap(:,3),val_EAC);
%             color = [c(1) c(2) c(3)];
%         else
%             val = trs_sig(7);
%             val_EAC = trs_EAC_sig(7);
%             color = [1 1 1];            
%         end
%         if n == 1
%             patch([0 0 1 1 ],[n_depth n_depth+1 n_depth+1 n_depth],color);
%         else
%             patch([n-1 n-1 n n],[n_depth n_depth+1 n_depth+1 n_depth],color);
%         end
%         x_pos = interp1(xlim_limits,xaxis_limits,n-0.6);
%         
%         if isnan(val)
%             val = 0;
%         end
%         if isnan(val_EAC)
%             val_EAC = 0;
%         end
%         
%         if sig_indicator == 1
%             annotation(gcf,'textbox',...
%             [x_pos y_pos_EAC 0.1 0.1],...
%             'String',[num2str(round(val,2)),' \color{black}{\diamondsuit}'],...
%             'LineStyle','none',...
%             'FontSize',14,...
%             'FitBoxToText','off');
%         else
%             annotation(gcf,'textbox',...
%             [x_pos y_pos_EAC 0.1 0.1],...
%             'String',[num2str(round(val,2))],...
%             'LineStyle','none',...
%             'FontSize',14,...
%             'FitBoxToText','off');                  
%         end
%         if sig_indicator_EAC == 1        
%             annotation(gcf,'textbox',...
%             [x_pos y_pos 0.1 0.1],...
%             'String',['\bf', num2str(round(val_EAC,2)),' \color{black}{\diamondsuit}'],...
%             'LineStyle','none',...
%             'FontSize',14,...
%             'FitBoxToText','off','FontName','Helvetica-Narrow');        
%         else
%             annotation(gcf,'textbox',...
%             [x_pos y_pos 0.1 0.1],...
%             'String',['\bf', num2str(round(val_EAC,2))],...
%             'LineStyle','none',...
%             'FontSize',14,...
%             'FitBoxToText','off','FontName','Helvetica-Narrow');     
%         end
%     end
% 
% end
% 
% cmap = cmocean('balance',23);
% cmap_limits = linspace(-1,1,23);
% xaxis_limits = linspace(0.06,0.9,6);
% yaxis_limits = linspace(0.11,0.93,6);
% xlim_limits = linspace(0,8,6);
% ylim_limits = linspace(1,8,6);
% 
% ylim([1 7])
% xlim([0 7])
% y_spacing = linspace(1.1,2,6);
% x_spacing = linspace(0,0.2,6);
% for n_depth = 1:6
% 
%     % organise trends over time
%     trs_sig = [trends_MAI(n_depth).t1950s_sig, trends_MAI(n_depth).t1960s_sig, trends_MAI(n_depth).t1970s_sig, trends_MAI(n_depth).t1980s_sig, ...
%         trends_MAI(n_depth).t1990s_sig, trends_MAI(n_depth).t2000s_sig, trends_MAI(n_depth).t2010s_sig, trends_MAI(n_depth).total_sig_1990s];
%     trs_insig = [trends_MAI(n_depth).t1950s_insig, trends_MAI(n_depth).t1960s_insig, trends_MAI(n_depth).t1970s_insig, trends_MAI(n_depth).t1980s_insig, ...
%         trends_MAI(n_depth).t1990s_insig, trends_MAI(n_depth).t2000s_insig, trends_MAI(n_depth).t2010s_insig, trends_MAI(n_depth).total_insig];
%     trs_EAC_sig = [trends_MAI(n_depth).t1950s_EAC_sig, trends_MAI(n_depth).t1960s_EAC_sig, trends_MAI(n_depth).t1970s_EAC_sig, trends_MAI(n_depth).t1980s_EAC_sig, ...
%         trends_MAI(n_depth).t1990s_EAC_sig, trends_MAI(n_depth).t2000s_EAC_sig, trends_MAI(n_depth).t2010s_EAC_sig, trends_MAI(n_depth).total_EAC_sig_1990s];    
%     trs_EAC_insig = [trends_MAI(n_depth).t1950s_EAC_insig, trends_MAI(n_depth).t1960s_EAC_insig, trends_MAI(n_depth).t1970s_EAC_insig, trends_MAI(n_depth).t1980s_EAC_insig, ...
%         trends_MAI(n_depth).t1990s_EAC_insig, trends_MAI(n_depth).t2000s_EAC_insig, trends_MAI(n_depth).t2010s_EAC_insig, trends_MAI(n_depth).total_EAC_insig];    
%     y_pos = interp1(ylim_limits,fliplr(yaxis_limits),n_depth+y_spacing(n_depth),'Linear','extrap');
%     y_pos_EAC = interp1(ylim_limits,fliplr(yaxis_limits),n_depth+y_spacing(n_depth)+0.3,'Linear','extrap');
%     xx = x_spacing(n_depth);
%     for n = 1:numel(trs_sig)-1
%         if n <= 7
%             % determine value used
%             if isfinite(trs_sig(n))
%                 val = trs_sig(n);
%                 sig_indicator = 1;
%             else
%                 val = trs_insig(n);
%                 sig_indicator = 0;
%             end
%             if isfinite(trs_EAC_sig(n))
%                 val_EAC = trs_EAC_sig(n);
%                 sig_indicator_EAC = 1;
%             else
%                 val_EAC = trs_EAC_insig(n);
%                 sig_indicator_EAC = 0;
%             end            
%             c(1) = interp1(cmap_limits,cmap(:,1),val_EAC);
%             c(2) = interp1(cmap_limits,cmap(:,2),val_EAC);
%             c(3) = interp1(cmap_limits,cmap(:,3),val_EAC);
%             color = [c(1) c(2) c(3)];
%             if sum(isfinite(color)) == 0
%                 color = [1 1 1];
%             end
%         else
%             val = trs_sig(7);
%             val_EAC = trs_EAC_sig(7);
%             color = [1 1 1];            
%         end
%         if n == 1
%             patch([0 0 1 1 ],[n_depth n_depth+1 n_depth+1 n_depth],color);
%         else
%             patch([n-1 n-1 n n],[n_depth n_depth+1 n_depth+1 n_depth],color);
%         end
%         x_pos = interp1(xlim_limits,xaxis_limits,n+xx);
%         
%         if isnan(val)
%             val = 0;
%         end
%         if isnan(val_EAC)
%             val_EAC = 0;
%         end
%         
%         if sig_indicator == 1
%             annotation(gcf,'textbox',...
%             [x_pos y_pos_EAC 0.1 0.1],...
%             'String',[num2str(round(val,2)),' \color{black}{\diamondsuit}'],...
%             'LineStyle','none',...
%             'FontSize',14,...
%             'FitBoxToText','off');
%         else
%             annotation(gcf,'textbox',...
%             [x_pos y_pos_EAC 0.1 0.1],...
%             'String',[num2str(round(val,2))],...
%             'LineStyle','none',...
%             'FontSize',14,...
%             'FitBoxToText','off');                  
%         end
%         if sig_indicator_EAC == 1        
%             annotation(gcf,'textbox',...
%             [x_pos y_pos 0.1 0.1],...
%             'String',['\bf', num2str(round(val_EAC,2)),' \color{black}{\diamondsuit}'],...
%             'LineStyle','none',...
%             'FontSize',14,...
%             'FitBoxToText','off','FontName','Helvetica-Narrow');        
%         else
%             annotation(gcf,'textbox',...
%             [x_pos y_pos 0.1 0.1],...
%             'String',['\bf', num2str(round(val_EAC,2))],...
%             'LineStyle','none',...
%             'FontSize',14,...
%             'FitBoxToText','off','FontName','Helvetica-Narrow');     
%         end
%     end
% 
% end
% 
% set(gca,'YDir','Reverse','LineWidth',1,'YTick',1.5:1:9.5,'YTickLabels',[2, 10, 20, 30, 40, 50, 85],...
%     'XTick',0.5:1:6.5,'XTickLabels',[{'1950s'} {'1960s'} {'1970s'} {'1980s'} {'1990s'} {'2000s'} {'2010s'}],'FontSize',16,'Box','On');
% ylabel('Depth [m]');
% title('NRS Maria Island')




