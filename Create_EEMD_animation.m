%% Load in data

addpath(genpath('C:\Users\mphem\Documents\Work\UNSW\'))
options.data_dir = 'C:\Users\mphem\Documents\Work\UNSW\Trends\Data\';
options.plot_dir = 'C:\Users\mphem\Documents\Work\UNSW\Trends\Plots\';

% NRSPHB
NRSPHB_data = load([options.data_dir,'NRSPHB_data.mat']);
NRSPHB_trends = load([options.data_dir,'NRSPHB_trends.mat']);

%% Convert time

% data time
nt = size(NRSPHB_data.t,2);
for nn = 1:numel(NRSPHB_trends.EEMD_trend)
    for t = 1:nt
        a = squeeze(vertcat(NRSPHB_data.t(nn,t,:)))';
        NRSPHB_data.t_conv(nn).t(t) = datenum(convertCharsToStrings(a));
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


%% Create animation

IMFs = NRSPHB_trends.EEMD_imfs.IMF_1;
IMF_t = NRSPHB_trends.EEMD_t_conv(n).t;

k=VideoWriter('C:\Users\mphem\Documents\Work\UNSW\Trends\Videos\EEMD_animation','MPEG-4');
open(k);

figure('units','normalized','position',[.1 .1 .6 .6]);
plot(IMF_t,sum(IMFs(:,:),1),'LineWidth',1.25)
set(gca,'LineWidth',3,'FontSize',16)
ylabel('Temperature Anomaly [^\circC]'); 
xlim([datenum(1950,01,01) datenum(2025,01,01)]);
ylim([-3.5 3.5]);
datetick('x','KeepLimits'); add_zero
%%%%%%%%%%%%%%
ax = axes
plot(IMF_t,sum(IMFs(:,:),1),'LineWidth',3,'Color','r')
set(gca,'visible','off');
xlim([datenum(1950,01,01) datenum(2025,01,01)]);
ylim([-3.5 3.5]);
F(1) = getframe;
%%%%%%%%%%%%%%
cla(ax)
ax = axes
plot(IMF_t,sum(IMFs(2:end,:),1),'LineWidth',3,'Color','r')
set(gca,'visible','off');
xlim([datenum(1950,01,01) datenum(2025,01,01)]);
ylim([-3.5 3.5]);
F(2) = getframe;
%%%%%%%%%%%%%%
cla(ax)
ax = axes
plot(IMF_t,sum(IMFs(3:end,:),1),'LineWidth',3,'Color','r')
set(gca,'visible','off');
xlim([datenum(1950,01,01) datenum(2025,01,01)]);
ylim([-3.5 3.5]);
F(3) = getframe;
%%%%%%%%%%%%%%
cla(ax)
ax = axes
plot(IMF_t,sum(IMFs(4:end,:),1),'LineWidth',3,'Color','r')
set(gca,'visible','off');
xlim([datenum(1950,01,01) datenum(2025,01,01)]);
ylim([-3.5 3.5]);
F(4) = getframe;
%%%%%%%%%%%%%%
cla(ax)
ax = axes
plot(IMF_t,sum(IMFs(5:end,:),1),'LineWidth',3,'Color','r')
set(gca,'visible','off');
xlim([datenum(1950,01,01) datenum(2025,01,01)]);
ylim([-3.5 3.5]);
F(5) = getframe;
%%%%%%%%%%%%%%
cla(ax)
ax = axes
plot(IMF_t,sum(IMFs(6:end,:),1),'LineWidth',3,'Color','r')
set(gca,'visible','off');
xlim([datenum(1950,01,01) datenum(2025,01,01)]);
ylim([-3.5 3.5]);
F(6) = getframe;
%%%%%%%%%%%%%%
cla(ax)
ax = axes
plot(IMF_t,sum(IMFs(7:end,:),1),'LineWidth',3,'Color','r')
set(gca,'visible','off');
xlim([datenum(1950,01,01) datenum(2025,01,01)]);
ylim([-3.5 3.5]);
F(7) = getframe;
%%%%%%%%%%%%%%
cla(ax)
ax = axes
plot(IMF_t,sum(IMFs(8:end,:),1),'LineWidth',3,'Color','r')
set(gca,'visible','off');
xlim([datenum(1950,01,01) datenum(2025,01,01)]);
ylim([-3.5 3.5]);
F(8) = getframe;
%%%%%%%%%%%%%%
cla(ax)
ax = axes
plot(IMF_t,sum(IMFs(7:end,:),1),'LineWidth',3,'Color','r')
set(gca,'visible','off');
xlim([datenum(1950,01,01) datenum(2025,01,01)]);
ylim([-3.5 3.5]);
F(9) = getframe;
%%%%%%%%%%%%%%
cla(ax)
ax = axes
plot(IMF_t,sum(IMFs(6:end,:),1),'LineWidth',3,'Color','r')
set(gca,'visible','off');
xlim([datenum(1950,01,01) datenum(2025,01,01)]);
ylim([-3.5 3.5]);
F(10) = getframe;
%%%%%%%%%%%%%%
cla(ax)
ax = axes
plot(IMF_t,sum(IMFs(5:end,:),1),'LineWidth',3,'Color','r')
set(gca,'visible','off');
xlim([datenum(1950,01,01) datenum(2025,01,01)]);
ylim([-3.5 3.5]);
F(11) = getframe;
%%%%%%%%%%%%%%
cla(ax)
ax = axes
plot(IMF_t,sum(IMFs(4:end,:),1),'LineWidth',3,'Color','r')
set(gca,'visible','off');
xlim([datenum(1950,01,01) datenum(2025,01,01)]);
ylim([-3.5 3.5]);
F(12) = getframe;
%%%%%%%%%%%%%%
cla(ax)
ax = axes
plot(IMF_t,sum(IMFs(3:end,:),1),'LineWidth',3,'Color','r')
set(gca,'visible','off');
xlim([datenum(1950,01,01) datenum(2025,01,01)]);
ylim([-3.5 3.5]);
F(13) = getframe;
%%%%%%%%%%%%%%
cla(ax)
ax = axes
plot(IMF_t,sum(IMFs(2:end,:),1),'LineWidth',3,'Color','r')
set(gca,'visible','off');
xlim([datenum(1950,01,01) datenum(2025,01,01)]);
ylim([-3.5 3.5]);
F(14) = getframe;
%%%%%%%%%%%%%%
cla(ax)
ax = axes
plot(IMF_t,sum(IMFs(1:end,:),1),'LineWidth',3,'Color','r')
set(gca,'visible','off');
xlim([datenum(1950,01,01) datenum(2025,01,01)]);
ylim([-3.5 3.5]);
F(15) = getframe;

writeVideo(k,F);
close(k);
close all