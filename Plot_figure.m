%% Load in data

addpath(genpath('C:\Users\mphem\Documents\Work\UNSW\'))
options.data_path = 'C:\Users\mphem\Documents\Work\UNSW\climatology\Revamped_scripts\Climatology\Data\';
options.plot_path = 'C:\Users\mphem\Documents\Work\UNSW\climatology\Revamped_scripts\Climatology\Plotting\';
data_path = 'C:\Users\mphem\Documents\Work\UNSW\climatology\Revamped_scripts\Climatology\Data\';
plot_path = 'C:\Users\mphem\Documents\Work\UNSW\Trends\Plots\';

load([data_path,'bathy_PH.mat']);
PHB.lon = x_grid;
PHB.lat = y_grid;
PHB.bathy = bathy_grid;

%% Create figure

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

