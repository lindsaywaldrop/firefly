%Figure edits
clear all
close all

% Some parameters for converting z scale
maxconc_air=0.0806;
maxconc_water=3.1282;
dx_air=0.0024;
dx_water=6.1035e-04;

% Concentration captured in water
openfig('ConcWater_2Dsurr_20180524.fig'); % Open surface fig from Yanyan
h=gcf;
axesObjs = get(h,'Children');
dataObjs = get(axesObjs,'Children');

water_x = dataObjs.XData;  % Extracts x data
water_y = dataObjs.YData;  % Extracts y data
water_z = dataObjs.ZData;  % Extracts z data

water_z = water_z.*(dx_water.^2)./maxconc_water;
close figure 1

openfig('ConAirVSConwater_2Dsurrogate.fig');
h2=gcf;
axesObjs2 = get(h2,'Children');
dataObjs2 = get(axesObjs2,'Children');

air_x = dataObjs2{1,1}.XData;  % Extracts x data
air_y = dataObjs2{1,1}.YData;  % Extracts y data
air_z = dataObjs2{1,1}.ZData;  % Extracts z data

air_z = air_z.*(dx_air.^2)./maxconc_air;

close figure 1



% Create figure
figure1 = figure('InvertHardcopy','off','Color',[1 1 1]);

% Create axes
axes1 = axes('Position',[0.13 0.11 0.334659090909091 0.815]);
hold(axes1,'on');

% Create surf
surf(water_x,water_y,water_z,'LineStyle','none');

% Create xlabel
xlabel({'Re'},'Margin',2);

% Create zlabel
zlabel({'Conc. Captured','in Water'},'Margin',2);

% Create ylabel
ylabel({'Angle'},'Margin',2);

% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[0 6]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[0 200]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim(axes1,[0 3]);
view(axes1,[-43.5 25.2]);
box(axes1,'on');
grid(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',16);
% Create axes
axes2 = axes('Position',[0.570340909090909 0.11 0.334659090909091 0.815]);
hold(axes2,'on');

% Create surf
surf(air_x,air_y,air_z,'LineStyle','none');

% Create xlabel
xlabel({'Gapwidth'},'Margin',2);

% Create zlabel
zlabel({'Conc. Captured ','in Air'},'Margin',2);

% Create ylabel
ylabel({'Angle'},'Margin',2);

% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes2,[0 60]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes2,[0 200]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim(axes2,[3000 7000]);
view(axes2,[-37.5 30]);
box(axes2,'on');
grid(axes2,'on');
% Set the remaining axes properties
set(axes2,'FontSize',16);

% Create textbox
annotation(figure1,'textbox',...
    [0.058153653599627 0.905441284490068 0.0419240953221536 0.0733590733590733],...
    'String',{'E'},...
    'LineStyle','none',...
    'FontSize',32,...
    'FontName','Arial');

% Create textbox
annotation(figure1,'textbox',...
    [0.499141043440525 0.913745797979333 0.0401588702559576 0.0733590733590733],...
    'String',{'F'},...
    'LineStyle','none',...
    'FontSize',32,...
    'FontName','Arial');