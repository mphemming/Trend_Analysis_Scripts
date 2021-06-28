%% Load in data

addpath(genpath('C:\Users\mphem\Documents\Work\UNSW\'))
options.data_dir = 'C:\Users\mphem\Documents\Work\UNSW\Trends\Data\';
options.plot_dir = 'C:\Users\mphem\Documents\Work\UNSW\Trends\Plots\';

% NRSPHB
NRSPHB_data = load([options.data_dir,'NRSPHB_data.mat']);
NRSPHB_trends = load([options.data_dir,'NRSPHB_trends.mat']);
NRSPHB_sig = load([options.data_dir,'NRSPHB_trends_server.mat']);
% NRSMAI
NRSMAI_data = load([options.data_dir,'NRSMAI_data.mat']);
NRSMAI_trends = load([options.data_dir,'NRSMAI_trends.mat']);
NRSMAI_sig = load(['C:\Users\mphem\Documents\Work\UNSW\Trends\Data\Split_data_server_deseason\NRSMAI_trends_server2.mat']);

%% Convert time

% data time
nt = size(NRSPHB_data.t,2);
for nn = 1:numel(NRSPHB_trends.EEMD_trend)
    for t = 1:nt
        a = squeeze(vertcat(NRSPHB_data.t(nn,t,:)))';
        NRSPHB_data.t_conv(nn).t(t) = datenum(convertCharsToStrings(a));
    end
end    
nt = size(NRSMAI_data.t,2);
for nn = 1:3
    for t = 1:nt
        a = squeeze(vertcat(NRSMAI_data.t(nn,t,:)))';
        NRSMAI_data.t_conv(nn).t(t) = datenum(convertCharsToStrings(a));
    end
end    

% EEMD time
for nn = 1:numel(NRSPHB_trends.EEMD_trend)
    a = NRSPHB_trends.EEMD_t{nn};
    for t = 1:size(a,1)
        b = a(t,:);
        NRSPHB_trends.EEMD_t_conv(nn).t(t) = datenum(convertCharsToStrings(b));
    end
end
for nn = 1:3
    a = squeeze(NRSMAI_trends.EEMD_t(nn,:,:));
    for t = 1:size(a,1)
        b = a(t,:);
        NRSMAI_trends.EEMD_t_conv(nn).t(t) = datenum(convertCharsToStrings(b));
    end
end

%% as MAI confidence starts in 1940 instead of 1944, need to shift

conf_t_min = nanmin(NRSMAI_data.t_conv(n).t);
trend_t_min = nanmin(NRSMAI_trends.EEMD_t_conv(n).t);
for n = 1:3
    NRSMAI_data.t_conv(n).t_shifted = NRSMAI_data.t_conv(n).t+ (trend_t_min-conf_t_min);
end


%% Create Monotonic trend plots

%% PHB

figure('units','normalized','position',[.1 .1 .6 .6]);
hold on;
for n = 1:7
    tr = NRSPHB_trends.EEMD_trend{n}; tr = tr-tr(1);
    % determine where significant
    clear sig
    for nn = 1:numel(tr)
        t = NRSPHB_trends.EEMD_t_conv(1).t(nn);
        f = find(NRSPHB_data.t_conv(1).t == t);
        if tr(nn) <= NRSPHB_sig.EEMD_conf_std_limit(1,f)
            sig(nn) = 0;
        else
            sig(nn) = 1;
        end
    end
    % plot
    plot(NRSPHB_trends.EEMD_t_conv(n).t,tr,'k','LineWidth',4)
    p(n) = plot(NRSPHB_trends.EEMD_t_conv(n).t,tr,'LineWidth',2)
    plot(NRSPHB_trends.EEMD_t_conv(n).t(sig == 0),tr(sig == 0),'--w','LineWidth',2)
end
leg = legend(p,[{'2'}, {'19'}, {'31'}, {'40'}, {'50'}, {'77'}, {'99'}])
set(leg,'FontSize',16,'Location','NorthWest','Box','Off');
set(gca,'FontSize',16,'LineWidth',2,'Box','On');
datetick
ylim([-0.05 0.8]); ylabel('Temperature Change [^\circ C]');
title('PHB - Monotonic trends');
grid on

print(gcf, '-dpng','-r400', ['C:\Users\mphem\Documents\Work\UNSW\Trends\Plots\','PHB_Monotonic_trends'])
close all

%% MAI

figure('units','normalized','position',[.1 .1 .6 .6]);
hold on;
clear p
for n = 1:3
    tr = NRSMAI_trends.EEMD_trend(n,:); tr = tr-tr(1);
    % determine where significant
    clear sig
    for nn = 1:numel(tr)
        t = NRSMAI_trends.EEMD_t_conv(1).t(nn);
        f = find(NRSMAI_data.t_conv(1).t == t);
        if tr(nn) <= NRSMAI_sig.EEMD_conf_std_limit(1,f)
            sig(nn) = 0;
        else
            sig(nn) = 1;
        end
    end
    % plot
    plot(NRSMAI_trends.EEMD_t_conv(n).t,tr,'k','LineWidth',4)
    p(n) = plot(NRSMAI_trends.EEMD_t_conv(n).t,tr,'LineWidth',2)
    plot(NRSMAI_trends.EEMD_t_conv(n).t(sig == 0),tr(sig == 0),'--w','LineWidth',2)
end
leg = legend(p,[{'2'}, {'20'}, {'50'}])
set(leg,'FontSize',16,'Location','NorthWest','Box','Off');
set(gca,'FontSize',16,'LineWidth',2,'Box','On');
xlim([datenum(1940,01,01) datenum(2025,01,01)])
datetick('x','KeepLimits');
ylim([-0.05 1]); ylabel('Temperature Change [^\circ C]');
title('MAI - Monotonic trends');
grid on

print(gcf, '-dpng','-r400', ['C:\Users\mphem\Documents\Work\UNSW\Trends\Plots\','MAI_Monotonic_trends'])
close all

%% Create Monotonic+IMF trend plots

%% PHB

figure('units','normalized','position',[.1 .1 .6 .6]);
hold on;
clear p
for n = 1:7
    tr = NRSPHB_trends.EEMD_trend_EAC{n}; tr = tr-tr(1);
    % determine where significant
    clear sig
    for nn = 1:numel(tr)
        t = NRSPHB_trends.EEMD_t_conv(1).t(nn);
        f = find(NRSPHB_data.t_conv(1).t == t);
        if tr(nn) <= NRSPHB_sig.EEMD_conf_std_limit_EAC(1,f)
            sig(nn) = 0;
        else
            sig(nn) = 1;
        end
    end
    % plot
    plot(NRSPHB_trends.EEMD_t_conv(n).t,tr,'k','LineWidth',4)
    p(n) = plot(NRSPHB_trends.EEMD_t_conv(n).t,tr,'LineWidth',2)
    plot(NRSPHB_trends.EEMD_t_conv(n).t(sig == 0),tr(sig == 0),'--w','LineWidth',2)
end
add_zero
leg = legend(p,[{'2'}, {'19'}, {'31'}, {'40'}, {'50'}, {'77'}, {'99'}])
set(leg,'FontSize',16,'Location','NorthWest','Box','Off');
set(gca,'FontSize',16,'LineWidth',2,'Box','On');
datetick
ylim([-0.5 1.3]); ylabel('Temperature Change [^\circ C]');
title('PHB - Monotonic+IMF_{last} trends');
grid on

print(gcf, '-dpng','-r400', ['C:\Users\mphem\Documents\Work\UNSW\Trends\Plots\','PHB_MonotonicIMF_trends'])
close all

%% MAI

figure('units','normalized','position',[.1 .1 .6 .6]);
hold on;
clear p
for n = 1:3
    tr = NRSMAI_trends.EEMD_trend_EAC(n,:); tr = tr-tr(1);
    % determine where significant
    clear sig
    for nn = 1:numel(tr)
        t = NRSMAI_trends.EEMD_t_conv(1).t(nn);
        f = find(NRSMAI_data.t_conv(1).t == t);
        if tr(nn) <= NRSMAI_sig.EEMD_conf_std_limit_EAC(1,f)
            sig(nn) = 0;
        else
            sig(nn) = 1;
        end
    end
    % plot
    plot(NRSMAI_trends.EEMD_t_conv(n).t,tr,'k','LineWidth',4)
    p(n) = plot(NRSMAI_trends.EEMD_t_conv(n).t,tr,'LineWidth',2)
    plot(NRSMAI_trends.EEMD_t_conv(n).t(sig == 0),tr(sig == 0),'--w','LineWidth',2)
end
leg = legend(p,[{'2'}, {'20'}, {'50'}])
set(leg,'FontSize',16,'Location','NorthWest','Box','Off');
set(gca,'FontSize',16,'LineWidth',2,'Box','On');
xlim([datenum(1940,01,01) datenum(2025,01,01)])
datetick
ylim([-0.05 1.6]); ylabel('Temperature Change [^\circ C]');
title('MAI - Monotonic+IMF_{last} trends');
grid on

print(gcf, '-dpng','-r400', ['C:\Users\mphem\Documents\Work\UNSW\Trends\Plots\','MAI_MonotonicIMF_trends'])
close all


