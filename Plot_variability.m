
%% load data

addpath(genpath('C:\Users\mphem\Documents\Work\UNSW\'))
main_path = '\Users\mphem\Documents\Work\UNSW\Trends\';
NRSPHB_agg = get_mooring([main_path,'Data\PH100_TEMP_1953-2020_aggregated_v1.nc'],1)
NRSMAI_agg = get_mooring([main_path,'Data\MAI090_TEMP_1944-2020_aggregated_v1.nc'],1)

%% get standard deviations

t_grid = datenum(1950:10:2020,01,01);
Ds = [2,20,50];

for n_t = 1:numel(t_grid)
    for n_D = 1:numel(Ds)
        %--------------------------------------------------------------
        % NRSPHB
        [yr,mn,dy] = datevec(NRSPHB_agg.TIME);
        c = NRSPHB_agg.TIME >= t_grid(n_t)-365.25*5 & NRSPHB_agg.TIME <= t_grid(n_t)+365.25*5 & ...
            NRSPHB_agg.DEPTH_AGG >= Ds(n_D)-2 & NRSPHB_agg.DEPTH_AGG <= Ds(n_D)+2;
        PHB_std(n_D,n_t) = nanstd(NRSPHB_agg.TEMP_AGG(c));
        PHB_n(n_D,n_t) = numel(NRSPHB_agg.TEMP_AGG(c)); 
        PHB_n_summer(n_D,n_t) = numel(NRSPHB_agg.TEMP_AGG(c & (mn == 12 | mn == 1 | mn == 2)));  % portion of point in summer
        PHB_n_autumn(n_D,n_t) = numel(NRSPHB_agg.TEMP_AGG(c & (mn == 3 | mn == 4 | mn == 5)));  % portion of point in autumn
        PHB_n_winter(n_D,n_t) = numel(NRSPHB_agg.TEMP_AGG(c & (mn == 6 | mn == 7 | mn == 8)));  % portion of point in winter
        PHB_n_spring(n_D,n_t) = numel(NRSPHB_agg.TEMP_AGG(c & (mn == 9 | mn == 10 | mn == 11)));  % portion of point in spring
        %--------------------------------------------------------------
        % NRSMAI
        [yr,mn,dy] = datevec(NRSMAI_agg.TIME);
        c = NRSMAI_agg.TIME >= t_grid(n_t)-365.25*5 & NRSMAI_agg.TIME <= t_grid(n_t)+365.25*5 & ...
            NRSMAI_agg.DEPTH_AGG >= Ds(n_D)-2 & NRSMAI_agg.DEPTH_AGG <= Ds(n_D)+2;
        MAI_std(n_D,n_t) = nanstd(NRSMAI_agg.TEMP_AGG(c));
        MAI_n(n_D,n_t) = numel(NRSMAI_agg.TEMP_AGG(c)); 
        MAI_n_summer(n_D,n_t) = numel(NRSMAI_agg.TEMP_AGG(c & (mn == 12 | mn == 1 | mn == 2)));  % portion of point in summer
        MAI_n_autumn(n_D,n_t) = numel(NRSMAI_agg.TEMP_AGG(c & (mn == 3 | mn == 4 | mn == 5)));  % portion of point in autumn
        MAI_n_winter(n_D,n_t) = numel(NRSMAI_agg.TEMP_AGG(c & (mn == 6 | mn == 7 | mn == 8)));  % portion of point in winter
        MAI_n_spring(n_D,n_t) = numel(NRSMAI_agg.TEMP_AGG(c & (mn == 9 | mn == 10 | mn == 11)));  % portion of point in spring        
    end
end

% need to create nice plot

%% binned profiles

t_grid = datenum(1950:10:2010,01,01);
D_PHB = [2, 19, 31, 40, 50, 77, 99];
D_MAI = [2, 20, 50];

for n_t = 1:numel(t_grid)
    % PHB
    c = NRSPHB_agg.TIME >= t_grid(n_t)-365.25*5 & NRSPHB_agg.TIME <= t_grid(n_t)+365.25*5;
    PHB_bins(n_t).profile = bin_variable_profile(NRSPHB_agg.TEMP_AGG(c),NRSPHB_agg.DEPTH_AGG(c),D_PHB,1);
    % MAI
    c = NRSMAI_agg.TIME >= t_grid(n_t)-365.25*5 & NRSMAI_agg.TIME <= t_grid(n_t)+365.25*5;
    MAI_bins(n_t).profile = bin_variable_profile(NRSMAI_agg.TEMP_AGG(c),NRSMAI_agg.DEPTH_AGG(c),D_MAI,1);    
end

%% NRSPHB figure
    
figure('units','normalized','position',[.1 .1 .7 .8]);  
 
% Create axes
axes('Parent',gcf,'Position',[0.13 0.11 0.386183035714286 0.815]);

cmap = cmocean('thermal',5);
nn = 0;
hold on
for n = [1:4,7]
    nn = nn+1;
    plot(PHB_bins(n).profile.mean_var_int,PHB_bins(n).profile.vertical,'LineWidth',4,'Color','k');
    p(n) = plot(PHB_bins(n).profile.mean_var_int,PHB_bins(n).profile.vertical,'LineWidth',2,'Color',cmap(nn,:));
    scatter(PHB_bins(n).profile.mean_var,PHB_bins(n).profile.vertical,'MarkerFaceColor',cmap(nn,:),'MarkerEdgeColor','k')
end

set(gca,'YDir','Reverse','LineWidth',2,'FontSize',16,'Box','On');
grid on;
leg = legend([p(1:4),p(7)],'1950s','1960s','1970s','1980','2010s')
set(leg,'Location','NorthWest','Box','Off')
title('Mean Profiles'); ylabel('Depth [m]'); xlabel('Temperature [^\circC]');
    
% Create axes
axes('Parent',gcf,...
    'Position',[0.537760416666667 0.11 0.367239583333333 0.815]);

cmap = cmocean('thermal',5);
nn = 0;
hold on
for n = [1:4,7]
    nn = nn+1;
    plot(PHB_bins(n).profile.mean_var_int-PHB_bins(1).profile.mean_var_int,PHB_bins(n).profile.vertical,'LineWidth',4,'Color','k');
    p(n) = plot(PHB_bins(n).profile.mean_var_int-PHB_bins(1).profile.mean_var_int,PHB_bins(n).profile.vertical,'LineWidth',2,'Color',cmap(nn,:));
    scatter(PHB_bins(n).profile.mean_var-PHB_bins(1).profile.mean_var_int,PHB_bins(n).profile.vertical,'MarkerFaceColor',cmap(nn,:),'MarkerEdgeColor','k')
end

set(gca,'YDir','Reverse','LineWidth',2,'FontSize',16,'Box','On','YTickLabels','');
grid on;
title('Temperature change since 1950s'); xlabel('Temperature [^\circC]');

print(gcf, '-dpng','-r400', ['C:\Users\mphem\Documents\Work\UNSW\Trends\Plots\','PHB_mean_profs'])
close all

%% NRSMAI figure
    
figure('units','normalized','position',[.1 .1 .7 .8]);  
 
% Create axes
axes('Parent',gcf,'Position',[0.13 0.11 0.386183035714286 0.815]);

cmap = cmocean('thermal',7);
nn = 0;
hold on
for n = 1:7
    nn = nn+1;
    plot(MAI_bins(n).profile.mean_var_int,MAI_bins(n).profile.vertical,'LineWidth',4,'Color','k');
    p(n) = plot(MAI_bins(n).profile.mean_var_int,MAI_bins(n).profile.vertical,'LineWidth',2,'Color',cmap(nn,:));
    scatter(MAI_bins(n).profile.mean_var,MAI_bins(n).profile.vertical,'MarkerFaceColor',cmap(nn,:),'MarkerEdgeColor','k')
end

set(gca,'YDir','Reverse','LineWidth',2,'FontSize',16,'Box','On','YLim',[0 100],'XLim',[12.5 15]);
grid on;
leg = legend(p,'1950s','1960s','1970s','1980','1990s','2000s','2010s')
set(leg,'Location','SouthWest','Box','Off')
title('Mean Profiles'); ylabel('Depth [m]'); xlabel('Temperature [^\circC]');


% Create axes
axes('Parent',gcf,...
    'Position',[0.537760416666667 0.11 0.367239583333333 0.815]);

cmap = cmocean('thermal',7);
nn = 0;
hold on
for n = 1:7
    nn = nn+1;
    plot(MAI_bins(n).profile.mean_var_int-MAI_bins(1).profile.mean_var_int,MAI_bins(n).profile.vertical,'LineWidth',4,'Color','k');
    p(n) = plot(MAI_bins(n).profile.mean_var_int-MAI_bins(1).profile.mean_var_int,MAI_bins(n).profile.vertical,'LineWidth',2,'Color',cmap(nn,:));
    scatter(MAI_bins(n).profile.mean_var-MAI_bins(1).profile.mean_var_int,MAI_bins(n).profile.vertical,'MarkerFaceColor',cmap(nn,:),'MarkerEdgeColor','k')
end

set(gca,'YDir','Reverse','LineWidth',2,'FontSize',16,'Box','On','YTickLabels','','YLim',[0 100]);
grid on;
title('Temperature change since 1950s'); xlabel('Temperature [^\circC]');

print(gcf, '-dpng','-r400', ['C:\Users\mphem\Documents\Work\UNSW\Trends\Plots\','MAI_mean_profs'])
close all

%% Mixed gradients

% PHB
% get unique dates
[yr,mn,dy] = datevec(NRSPHB_agg.TIME);
PHB_results.dates = datenum(yr,mn,dy);
PHB_results.un_dates = unique(PHB_results.dates);
mooring_first_date = nanmin(PHB_results.dates(NRSPHB_agg.TEMP_DATA_PLATFORM_AGG == 3 | NRSPHB_agg.TEMP_DATA_PLATFORM_AGG == 4));
% get dates where satellite for randomising profiles below
PHB_results.sat_dates = unique(PHB_results.dates(NRSPHB_agg.TEMP_DATA_PLATFORM_AGG == 4));
% randomise mooring dates for similar resolution with bottles and CTDs
PHB_results.un_dates_mooring = PHB_results.un_dates(PHB_results.un_dates >= mooring_first_date);
PHB_results.un_dates = [PHB_results.un_dates(PHB_results.un_dates < mooring_first_date)' PHB_results.sat_dates(1:3:end)'];
% get gradient between surface and 20m for each unique date
for n_date = 1:numel(PHB_results.un_dates)
    c_surf = PHB_results.dates == PHB_results.un_dates(n_date) & NRSPHB_agg.DEPTH_AGG < 4;
    c_31 = PHB_results.dates == PHB_results.un_dates(n_date) & NRSPHB_agg.DEPTH_AGG >= 29 & NRSPHB_agg.DEPTH_AGG < 33;
    if sum(c_surf) ~= 0 & sum(c_31) ~= 0
        PHB_results.Tdiff(n_date) = abs(nanmedian(NRSPHB_agg.TEMP_AGG(c_surf))-nanmedian(NRSPHB_agg.TEMP_AGG(c_31)));
        PHB_results.DOY(n_date) = datenum(0,unique(mn(c_surf)),unique(dy(c_31)));
        if PHB_results.Tdiff(n_date) <0.2 
%             if PHB_results.Tdiff(n_date) == 0 
%                 disp(n_date)
%             end
            PHB_results.mixed(n_date) = 1;
        else
            PHB_results.mixed(n_date) = 0;
        end
    else
        PHB_results.Tdiff(n_date) = NaN;
        PHB_results.DOY(n_date) = NaN;
        PHB_results.mixed(n_date) = NaN;
    end
end
% for each year proportion of data mixed
% NRSPHB
[yr,mn,dy] = datevec(PHB_results.un_dates);
n = 0;
for yr_n= 1950:10:2010
    n=n+1;
    c = yr >= yr_n & yr < yr_n+10 & isfinite(PHB_results.Tdiff) & PHB_results.Tdiff ~= 0;
    PHB_results.prop_mixed(n) = numel(PHB_results.Tdiff(c & PHB_results.mixed == 1))/numel(PHB_results.Tdiff(c))*100;
    PHB_results.prop_mixed_summer(n) = numel(PHB_results.Tdiff(c & PHB_results.mixed == 1  & ( mn == 12 | mn == 1 | mn == 2)))/numel(PHB_results.Tdiff(c))*100;
    PHB_results.prop_mixed_autumn(n) = numel(PHB_results.Tdiff(c & PHB_results.mixed == 1  & ( mn == 3 | mn == 4 | mn == 5)))/numel(PHB_results.Tdiff(c))*100;
    PHB_results.prop_mixed_winter(n) = numel(PHB_results.Tdiff(c & PHB_results.mixed == 1  & ( mn == 6 | mn == 7 | mn == 8)))/numel(PHB_results.Tdiff(c))*100;
    PHB_results.prop_mixed_spring(n) = numel(PHB_results.Tdiff(c & PHB_results.mixed == 1  & ( mn == 9 | mn == 10 | mn == 11)))/numel(PHB_results.Tdiff(c))*100;
    PHB_results.prop_n(n) = numel(PHB_results.Tdiff(c & isfinite(PHB_results.Tdiff)));
    PHB_results.prop_n_summer(n) = numel(PHB_results.Tdiff(c & ( mn == 12 | mn == 1 | mn == 2)));
    PHB_results.prop_n_autumn(n) = numel(PHB_results.Tdiff(c & ( mn == 3 | mn == 4 | mn == 5)));    
    PHB_results.prop_n_winter(n) = numel(PHB_results.Tdiff(c & ( mn == 6 | mn == 7 | mn == 8)));    
    PHB_results.prop_n_spring(n) = numel(PHB_results.Tdiff(c & ( mn == 9 | mn == 10 | mn == 11)));    
    PHB_results.mean_Tdiff(n) = nanmean(PHB_results.Tdiff(c));
    PHB_results.std_Tdiff(n) = nanstd(PHB_results.Tdiff(c));
    
    % some QC (if biased towards one season, if not enough data)
    % zeros in Tdiff seem legit (checked)
    if (PHB_results.prop_n_summer(n)/PHB_results.prop_n(n) <0.1 | ...
            PHB_results.prop_n_winter(n)/PHB_results.prop_n(n) <0.1 | ...
                PHB_results.prop_n_autumn(n)/PHB_results.prop_n(n) <0.1 | ...
                    PHB_results.prop_n_spring(n)/PHB_results.prop_n(n) <0.1) | ...
                        PHB_results.prop_n(n) < 150
        PHB_results.prop_mixed(n) = NaN;
        PHB_results.prop_mixed_summer(n) = NaN;
        PHB_results.prop_mixed_autumn(n) = NaN;
        PHB_results.prop_mixed_winter(n) = NaN;
        PHB_results.prop_mixed_spring(n) = NaN;
        PHB_results.mean_Tdiff(n) = NaN;
        PHB_results.std_Tdiff(n) = NaN;
    end

end



% MAI
% get unique dates
[yr,mn,dy] = datevec(NRSMAI_agg.TIME);
MAI_results.dates = datenum(yr,mn,dy);
MAI_results.un_dates = unique(MAI_results.dates);
mooring_first_date = nanmin(MAI_results.dates(NRSMAI_agg.TEMP_DATA_PLATFORM_AGG == 3 | NRSMAI_agg.TEMP_DATA_PLATFORM_AGG == 4));
% get dates where satellite for randomising profiles below
MAI_results.sat_dates = unique(MAI_results.dates(NRSMAI_agg.TEMP_DATA_PLATFORM_AGG == 4));
% randomise mooring dates for similar resolution with bottles and CTDs
MAI_results.un_dates_mooring = MAI_results.un_dates(MAI_results.un_dates >= mooring_first_date);
MAI_results.un_dates = [MAI_results.un_dates(MAI_results.un_dates < mooring_first_date)' MAI_results.sat_dates(1:1:end)'];
% get gradient between surface and 20m for each unique date
for n_date = 1:numel(MAI_results.un_dates)
    c_surf = MAI_results.dates == MAI_results.un_dates(n_date) & NRSMAI_agg.DEPTH_AGG < 4;
    c_31 = MAI_results.dates == MAI_results.un_dates(n_date) & NRSMAI_agg.DEPTH_AGG >= 29 & NRSMAI_agg.DEPTH_AGG < 33;
    if sum(c_surf) ~= 0 & sum(c_31) ~= 0
        MAI_results.Tdiff(n_date) = abs(nanmedian(NRSMAI_agg.TEMP_AGG(c_surf))-nanmedian(NRSMAI_agg.TEMP_AGG(c_31)));
        MAI_results.DOY(n_date) = datenum(0,unique(mn(c_surf)),unique(dy(c_31)));
        if MAI_results.Tdiff(n_date) <0.2 
%             if MAI_results.Tdiff(n_date) == 0 
%                 disp(n_date)
%             end
            MAI_results.mixed(n_date) = 1;
        else
            MAI_results.mixed(n_date) = 0;
        end
    else
        MAI_results.Tdiff(n_date) = NaN;
        MAI_results.DOY(n_date) = NaN;
        MAI_results.mixed(n_date) = NaN;
    end
end
% for each year proportion of data mixed
% NRSMAI
[yr,mn,dy] = datevec(MAI_results.un_dates);
n = 0;
for yr_n= 1950:10:2010
    n=n+1;
    c = yr >= yr_n & yr < yr_n+10 & isfinite(MAI_results.Tdiff) & MAI_results.Tdiff ~= 0;
    MAI_results.prop_mixed(n) = numel(MAI_results.Tdiff(c & MAI_results.mixed == 1))/numel(MAI_results.Tdiff(c))*100;
    MAI_results.prop_mixed_summer(n) = numel(MAI_results.Tdiff(c & MAI_results.mixed == 1  & ( mn == 12 | mn == 1 | mn == 2)))/numel(MAI_results.Tdiff(c))*100;
    MAI_results.prop_mixed_autumn(n) = numel(MAI_results.Tdiff(c & MAI_results.mixed == 1  & ( mn == 3 | mn == 4 | mn == 5)))/numel(MAI_results.Tdiff(c))*100;
    MAI_results.prop_mixed_winter(n) = numel(MAI_results.Tdiff(c & MAI_results.mixed == 1  & ( mn == 6 | mn == 7 | mn == 8)))/numel(MAI_results.Tdiff(c))*100;
    MAI_results.prop_mixed_spring(n) = numel(MAI_results.Tdiff(c & MAI_results.mixed == 1  & ( mn == 9 | mn == 10 | mn == 11)))/numel(MAI_results.Tdiff(c))*100;
    MAI_results.prop_n(n) = numel(MAI_results.Tdiff(c & isfinite(MAI_results.Tdiff)));
    MAI_results.prop_n_summer(n) = numel(MAI_results.Tdiff(c & ( mn == 12 | mn == 1 | mn == 2)));
    MAI_results.prop_n_autumn(n) = numel(MAI_results.Tdiff(c & ( mn == 3 | mn == 4 | mn == 5)));    
    MAI_results.prop_n_winter(n) = numel(MAI_results.Tdiff(c & ( mn == 6 | mn == 7 | mn == 8)));    
    MAI_results.prop_n_spring(n) = numel(MAI_results.Tdiff(c & ( mn == 9 | mn == 10 | mn == 11)));    
    MAI_results.mean_Tdiff(n) = nanmean(MAI_results.Tdiff(c));
    MAI_results.std_Tdiff(n) = nanstd(MAI_results.Tdiff(c));
    
    % some QC (if biased towards one season, if not enough data)
    % zeros in Tdiff seem legit (checked)
    if (MAI_results.prop_n_summer(n)/MAI_results.prop_n(n) <0.1 | ...
            MAI_results.prop_n_winter(n)/MAI_results.prop_n(n) <0.1 | ...
                MAI_results.prop_n_autumn(n)/MAI_results.prop_n(n) <0.1 | ...
                    MAI_results.prop_n_spring(n)/MAI_results.prop_n(n) <0.1)
        MAI_results.prop_mixed(n) = NaN;
        MAI_results.prop_mixed_summer(n) = NaN;
        MAI_results.prop_mixed_autumn(n) = NaN;
        MAI_results.prop_mixed_winter(n) = NaN;
        MAI_results.prop_mixed_spring(n) = NaN;
        MAI_results.mean_Tdiff(n) = NaN;
        MAI_results.std_Tdiff(n) = NaN;
    end

end



%% Figure mean change and std temperature difference between 1950s and 2010s

%% NRSPHB
figure('units','normalized','position',[.1 .1 .6 .6]);

axes('Parent',gcf,...
    'Position',[0.13 0.405478395061728 0.775 0.519521604938266]);
p1 = plot(1:7,inpaint_nans(PHB_results.mean_Tdiff,4), 'LineWidth',2); hold on;
scatter(1:7,PHB_results.mean_Tdiff, 'MarkerFaceColor',[0 0.447 0.741], 'MarkerEdgeColor',[0 0.447 0.741]);
p2 = plot(1:7,inpaint_nans(PHB_results.std_Tdiff,4), 'LineWidth',2);
scatter(1:7,PHB_results.std_Tdiff, 'MarkerFaceColor',[.929 .694 .125], 'MarkerEdgeColor',[.929 .694 .125]);
set(gca,'XLim',[0 8],'XTickLabels','','FontSize',16,'Box','On','LineWidth',2,'XGrid','On')
ylabel('Temperature [^\circC]')

axes('Parent',gcf,...
    'Position',[0.13 0.405478395061728 0.775 0.519521604938266]);
hold on
p3 = plot(1:7,inpaint_nans(PHB_results.prop_mixed,4), 'LineWidth',2, 'Color',[.466 .674 .188]);
scatter(1:7,PHB_results.prop_mixed, 'MarkerFaceColor',[.466 .674 .188], 'MarkerEdgeColor',[.466 .674 .188]);
set(gca,'YLim',[0 60],'Visible','Off','XLim',[0 8]);

axes('Parent',gcf,...
    'Position',[0.962456597222222 0.405478395061727 0.00217013888888884 0.519521604938264]);
set(gca,'YLim',[0 60],'Color',[.466 .674 .188],'FontSize',16,'YColor',[.466 .674 .188],'LineWidth',2);

leg = legend([p1 p2 p3],'Mean|T_{surf}  -  T_{31m}|','SD|T_{surf}  -  T_{31m}|','% profiles mixed')
set(leg,'Location','NorthWest','Box','Off');


axes('Parent',gcf,...
    'Position',[0.13 0.0860339506172837 0.775 0.308641975308642]);

y(1,:) = PHB_results.prop_n_summer./PHB_results.prop_n*100;
y(2,:) = PHB_results.prop_n_autumn./PHB_results.prop_n*100;
y(3,:) = PHB_results.prop_n_winter./PHB_results.prop_n*100;
y(4,:) = PHB_results.prop_n_spring./PHB_results.prop_n*100;

b = bar(y','stacked')
set(gca,'Xtick',1:7,'XTickLabels',[{'1950s'}, {'1960s'}, {'1970s'}, {'1980s'}, {'1990s'}, {'2000s'}, {'2010s'}],...
    'FontSize',16,'Box','On','LineWidth',2,'XLim',[0 8],'YTickLabels','')
set(b(1),'FaceColor',[.9 .7 0]);
set(b(2),'FaceColor',[.6 .4 .4]);
set(b(3),'FaceColor',[.4 .6 .6]);
set(b(4),'FaceColor',[0 .8 .5]);
leg = legend([b(1) b(2) b(3) b(4)],'Sum','Aut','Win','Spr')
set(leg,'Position',[0.0598234943355679 0.137704512230447 0.0987413204792473 0.201967597154923]);

axes('Parent',gcf,...
    'Position',[0.13 0.0860339506172837 0.775 0.308641975308642]);
plot(PHB_results.prop_n,'--k','LineWidth',2)
set(gca,'XLim',[0 8],'Visible','Off','YLim',[nanmin(PHB_results.prop_n) nanmax(PHB_results.prop_n)]);

axes1 = axes('Parent',gcf,...
    'Position',[0.962491319444445 0.0860339506172837 0.00217013888888884 0.308641975308642]);
set(gca,'YLim',[nanmin(PHB_results.prop_n) nanmax(PHB_results.prop_n)],'FontSize',16,'Box','On','LineWidth',2);

print(gcf, '-dpng','-r400', ['C:\Users\mphem\Documents\Work\UNSW\Trends\Plots\','PHB_prop_mean_std_change_time'])
close all

%% NRSMAI
figure('units','normalized','position',[.1 .1 .6 .6]);

axes('Parent',gcf,...
    'Position',[0.13 0.405478395061728 0.775 0.519521604938266]);
p1 = plot(1:7,inpaint_nans(MAI_results.mean_Tdiff,4), 'LineWidth',2); hold on;
scatter(1:7,MAI_results.mean_Tdiff, 'MarkerFaceColor',[0 0.447 0.741], 'MarkerEdgeColor',[0 0.447 0.741]);
p2 = plot(1:7,inpaint_nans(MAI_results.std_Tdiff,4), 'LineWidth',2);
scatter(1:7,MAI_results.std_Tdiff, 'MarkerFaceColor',[.929 .694 .125], 'MarkerEdgeColor',[.929 .694 .125]);
set(gca,'XLim',[0 8],'XTickLabels','','FontSize',16,'Box','On','LineWidth',2,'XGrid','On','YLim',[0 1.1])
ylabel('Temperature [^\circC]')

axes('Parent',gcf,...
    'Position',[0.13 0.405478395061728 0.775 0.519521604938266]);
hold on
p3 = plot(1:7,inpaint_nans(MAI_results.prop_mixed,4), 'LineWidth',2, 'Color',[.466 .674 .188]);
scatter(1:7,MAI_results.prop_mixed, 'MarkerFaceColor',[.466 .674 .188], 'MarkerEdgeColor',[.466 .674 .188]);
set(gca,'YLim',[30 80],'Visible','Off','XLim',[0 8]);

axes('Parent',gcf,...
    'Position',[0.962456597222222 0.405478395061727 0.00217013888888884 0.519521604938264]);
set(gca,'YLim',[30 80],'Color',[.466 .674 .188],'FontSize',16,'YColor',[.466 .674 .188],'LineWidth',2);

leg = legend([p1 p2 p3],'Mean|T_{surf}  -  T_{31m}|','SD|T_{surf}  -  T_{31m}|','% profiles mixed')
set(leg,'Location','NorthWest','Box','Off');


axes('Parent',gcf,...
    'Position',[0.13 0.0860339506172837 0.775 0.308641975308642]);

y(1,:) = MAI_results.prop_n_summer./MAI_results.prop_n*100;
y(2,:) = MAI_results.prop_n_autumn./MAI_results.prop_n*100;
y(3,:) = MAI_results.prop_n_winter./MAI_results.prop_n*100;
y(4,:) = MAI_results.prop_n_spring./MAI_results.prop_n*100;

b = bar(y','stacked')
set(gca,'Xtick',1:7,'XTickLabels',[{'1950s'}, {'1960s'}, {'1970s'}, {'1980s'}, {'1990s'}, {'2000s'}, {'2010s'}],...
    'FontSize',16,'Box','On','LineWidth',2,'XLim',[0 8],'YTickLabels','')
set(b(1),'FaceColor',[.9 .7 0]);
set(b(2),'FaceColor',[.6 .4 .4]);
set(b(3),'FaceColor',[.4 .6 .6]);
set(b(4),'FaceColor',[0 .8 .5]);
leg = legend([b(1) b(2) b(3) b(4)],'Sum','Aut','Win','Spr')
set(leg,'Position',[0.0598234943355679 0.137704512230447 0.0987413204792473 0.201967597154923]);

axes('Parent',gcf,...
    'Position',[0.13 0.0860339506172837 0.775 0.308641975308642]);
plot(MAI_results.prop_n,'--k','LineWidth',2)
set(gca,'XLim',[0 8],'Visible','Off','YLim',[nanmin(MAI_results.prop_n)-10 nanmax(MAI_results.prop_n)+10]);

axes1 = axes('Parent',gcf,...
    'Position',[0.962491319444445 0.0860339506172837 0.00217013888888884 0.308641975308642]);
set(gca,'YLim',[nanmin(MAI_results.prop_n)-10 nanmax(MAI_results.prop_n)+10],'FontSize',16,'Box','On','LineWidth',2);

print(gcf, '-dpng','-r400', ['C:\Users\mphem\Documents\Work\UNSW\Trends\Plots\','MAI_prop_mean_std_change_time'])
close all


