%% Load in data

addpath(genpath('C:\Users\mphem\Documents\Work\UNSW\'))
options.data_path = 'C:\Users\mphem\Documents\Work\UNSW\climatology\Revamped_scripts\Climatology\Data\';
options.plot_path = 'C:\Users\mphem\Documents\Work\UNSW\climatology\Revamped_scripts\Climatology\Plotting\';
data_path = 'C:\Users\mphem\Documents\Work\UNSW\climatology\Revamped_scripts\Climatology\Data\';
plot_path = 'C:\Users\mphem\Documents\Work\UNSW\Trends\Plots\';

load(['C:\Users\mphem\Documents\Work\UNSW\BATHYMETRY\','bathy_eac_etopo1.mat']);
% PHB.lon = x_grid;
% PHB.lat = y_grid;
% PHB.bathy = bathy_grid;

% SSTAARS climatology
filename = 'C:\Users\mphem\Documents\Work\UNSW\Trends\Data\SSTAARS.nc'

TEMP_mean = ncread(filename,'TEMP_mean');  
TEMP_mean_std_err = ncread(filename,'TEMP_mean_std_err');  
TEMP_trend = ncread(filename,'TEMP_trend');  
% get trend per decade
TEMP_trend = (TEMP_trend / numel(datenum(1992,03,22):1:datenum(2016,12,31)))*(365*10);

TEMP_trend_std_err = ncread(filename,'TEMP_trend_std_err');  
LONGITUDE = ncread(filename,'LONGITUDE');    
LATITUDE = ncread(filename,'LATITUDE');    
% LONG_mat = repmat(LONGITUDE,[1 32]);
% LAT_mat = repmat(LATITUDE,[1 92])';
[X,Y] = meshgrid(LONGITUDE,LATITUDE);


% get interpolated higher resolution
% [X Y] = meshgrid(nanmin(LONGITUDE):0.05:nanmax(LONGITUDE),nanmin(LATITUDE):0.05:nanmax(LATITUDE));
% Tm = TEMP_mean(:,:,1);
% T = griddata(double(LONG_mat(:)),double(LAT_mat(:)),double(Tm(:)),double(X),double(Y));


API = load('C:\Users\mphem\Documents\Work\UNSW\climatology\Revamped_scripts\Climatology\Utilities\Plot_google_map\api_key');
cd 'C:\Users\mphem\Documents\Work\UNSW\climatology\Revamped_scripts\Climatology\Utilities\'
  
% the reference stations
Locations.NRSPHB = [151.216 -34.118];
Locations.NRSMAI = [148.23 -42.6];
Locations.NRSNSI = [153.562 -27.342]
Locations.BMP120 = [150.3134 -36.2064];
Locations.CH100 = [153.3957 -30.2656];

%% get land mask using google earth
% xlim([111 157])
% ylim([-45 -15])
% [lonVec latVec imag] = plot_google_map('APIKey',API.apiKey,'ShowLabels',0);
% mask = ones(size(imag(:,:,1)));
% mask(imag(:,:,1) == 170) = 1;
% mask(imag(:,:,1) ~= 170) = 0;
% lv = repmat(lonVec,[1280,1]);
% lav = repmat(latVec,[1280,1])';
% m = griddata(lv(:),lav(:),mask(:),double(X),double(Y));
% 
% T(m == 0) = NaN;

%% Create figure

figure('units','normalized','position',[.1 0 .45 .9])
xlim([146 156])
ylim([-44 -25])
plot_google_map('APIKey',API.apiKey,'MapType','satellite','Resize',2,'ScaleWidth',0.2,'FigureResizeUpdate',0)
hold on
contourf(X,Y,TEMP_trend',100,'LineStyle','None')
xlim([146 156])
ylim([-44 -25])

cm = cbrewer('div', 'Spectral',30);
colormap(flipud(cm));
% caxis([12 24]);
caxis([0 0.5]);

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

cb = colorbar;
ylabel(cb,'Surface Temperature Trend [\circC decade^{-1}]');
xlabel('Longitude [^\circ E]')
ylabel('Latitude [^\circ S]')
set(gca,'FontSize',20,'LineWidth',2,'Box','On','YTick',[-42 -38 -34 -30 -26],'YTickLabels',[{'42'} {'38'} {'34'} {'30'} {'26'}]);

print(gcf, '-dpng','-r400', [plot_path,'Map_1'])

% figure('units','normalized','position',[.1 0 .6 .7])

% style = ['https://maps.googleapis.com/maps/api/staticmap?key=YOUR_API_KEY&center=-29.255517051688706,133.41701195000005&zoom=4&format=png&maptype=roadmap&style=element:geometry%7Ccolor:0xf5f5f5&style=element:geometry.fill%7Ccolor:0xffffff&style=element:labels%7Cvisibility:off&style=element:labels.icon%7Cvisibility:off&style=element:labels.text.fill%7Ccolor:0x616161&style=element:labels.text.stroke%7Ccolor:0xf5f5f5&style=feature:administrative%7Celement:geometry%7Cvisibility:off&style=feature:administrative.land_parcel%7Cvisibility:off&style=feature:administrative.land_parcel%7Celement:labels.text.fill%7Ccolor:0xbdbdbd&style=feature:administrative.neighborhood%7Cvisibility:off&style=feature:poi%7Cvisibility:off&style=feature:poi%7Celement:geometry%7Ccolor:0xeeeeee&style=feature:poi%7Celement:labels.text.fill%7Ccolor:0x757575&style=feature:poi.park%7Celement:geometry%7Ccolor:0xe5e5e5&style=feature:poi.park%7Celement:labels.text.fill%7Ccolor:0x9e9e9e&style=feature:road%7Cvisibility:off&style=feature:road%7Celement:geometry%7Ccolor:0xffffff&style=feature:road%7Celement:labels.icon%7Cvisibility:off&style=feature:road.arterial%7Celement:labels.text.fill%7Ccolor:0x757575&style=feature:road.highway%7Celement:geometry%7Ccolor:0xdadada&style=feature:road.highway%7Celement:labels.text.fill%7Ccolor:0x616161&style=feature:road.local%7Celement:labels.text.fill%7Ccolor:0x9e9e9e&style=feature:transit%7Cvisibility:off&style=feature:transit.line%7Celement:geometry%7Ccolor:0xe5e5e5&style=feature:transit.station%7Celement:geometry%7Ccolor:0xeeeeee&style=feature:water%7Celement:geometry%7Ccolor:0xc9c9c9&style=feature:water%7Celement:labels.text.fill%7Ccolor:0x9e9e9e&size=480x360']
% 
% xlim([111 155]);
% ylim([-45 -11]);
% plot_google_map('Style',style);
% hold on;
% 
% print(gcf, '-dpng','-r400', [plot_path,'Map_2'])













API = load('C:\Users\mphem\Documents\Work\UNSW\climatology\Revamped_scripts\Climatology\Utilities\Plot_google_map\api_key');

style = ['https://maps.googleapis.com/maps/api/staticmap?key=YOUR_API_KEY&center=-33.95478802426896,151.1297292205634&zoom=11&format=png&maptype=roadmap&style=element:geometry%7Ccolor:0x1d2c4d&style=element:labels%7Cvisibility:off&style=element:labels.text.fill%7Ccolor:0x8ec3b9&style=element:labels.text.stroke%7Ccolor:0x1a3646&style=feature:administrative%7Celement:geometry%7Cvisibility:off&style=feature:administrative.country%7Celement:geometry.stroke%7Ccolor:0x4b6878&style=feature:administrative.land_parcel%7Cvisibility:off&style=feature:administrative.land_parcel%7Celement:labels.text.fill%7Ccolor:0x64779e&style=feature:administrative.neighborhood%7Cvisibility:off&style=feature:administrative.province%7Celement:geometry.stroke%7Ccolor:0x4b6878&style=feature:landscape%7Celement:geometry.stroke%7Cweight:0.5&style=feature:landscape%7Celement:labels.icon%7Clightness:5&style=feature:landscape.man_made%7Ccolor:0x667d8f%7Csaturation:-10&style=feature:landscape.man_made%7Celement:geometry.stroke%7Ccolor:0x334e87&style=feature:landscape.natural%7Celement:geometry%7Ccolor:0x667d8f&style=feature:poi%7Cvisibility:off&style=feature:poi%7Celement:geometry%7Ccolor:0x283d6a&style=feature:poi%7Celement:labels.text.fill%7Ccolor:0x6f9ba5&style=feature:poi%7Celement:labels.text.stroke%7Ccolor:0x1d2c4d&style=feature:poi.park%7Celement:geometry.fill%7Ccolor:0x023e58&style=feature:poi.park%7Celement:labels.text.fill%7Ccolor:0x3C7680&style=feature:road%7Cvisibility:off&style=feature:road%7Celement:geometry%7Ccolor:0x304a7d&style=feature:road%7Celement:labels.icon%7Cvisibility:off&style=feature:road%7Celement:labels.text.fill%7Ccolor:0x98a5be&style=feature:road%7Celement:labels.text.stroke%7Ccolor:0x1d2c4d&style=feature:road.highway%7Celement:geometry%7Ccolor:0x2c6675&style=feature:road.highway%7Celement:geometry.stroke%7Ccolor:0x255763&style=feature:road.highway%7Celement:labels.text.fill%7Ccolor:0xb0d5ce&style=feature:road.highway%7Celement:labels.text.stroke%7Ccolor:0x023e58&style=feature:transit%7Cvisibility:off&style=feature:transit%7Celement:labels.text.fill%7Ccolor:0x98a5be&style=feature:transit%7Celement:labels.text.stroke%7Ccolor:0x1d2c4d&style=feature:transit.line%7Celement:geometry.fill%7Ccolor:0x283d6a&style=feature:transit.station%7Celement:geometry%7Ccolor:0x3a4762&style=feature:water%7Ccolor:0xe7e9ee&style=feature:water%7Celement:geometry%7Ccolor:0xe9edf7&style=feature:water%7Celement:labels.text.fill%7Ccolor:0x4e6d70&size=480x360'];

figure('units','normalized','position',[0 0.1 .5 .8]);
xlim([150.9 151.4]);
ylim([-34.3 -33.6]);
plot_google_map('Style',style);
hold on;

scatter(151.216, -34.118,150,'MarkerFaceColor','k','MarkerEdgeColor','k','Marker','Sq');
s4 = scatter(151.216, -34.118,60,'MarkerFaceColor','w','MarkerEdgeColor','k','Marker','Sq');

set(gca,'FontSize',22,'LineWidth',2,'Box','On');
xlabel('Longitude [^\circ E]');
ylabel('Latitude [^\circ S]');

print(gcf,'-dpng','-r400',[plot_path,'PHB_map']);

%% Australia-wide figure


figure('units','normalized','position',[0 0.1 .5 .5]);
xlim([109 159]);
ylim([-45 -5]);
plot_google_map('MapType','satellite');
hold on;

scatter(151.216, -34.118,150,'MarkerFaceColor','k','MarkerEdgeColor','k','Marker','Sq');

set(gca,'YTickLabels','','XTickLabels','','LineWidth',5,'Box','On');

print(gcf,'-dpng','-r400',[plot_path,'Australia_PHB']);

scatter(148.23, -42.6,150,'MarkerFaceColor','k','MarkerEdgeColor','k','Marker','Sq');
s4 = scatter(148.23, -42.6,60,'MarkerFaceColor','w','MarkerEdgeColor','k','Marker','Sq');

