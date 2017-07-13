function plot_concentration_v2(ind)
% ind =8;

%% Current file directory
file_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', '..');

%% dlsode time
filename = fullfile(file_dir, 'output', 'time_dlsode_fraction.csv');
delimiter = '';
formatSpec = '%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
fclose(fileID);
dlsode_time = dataArray{:, 1};
clearvars filename delimiter formatSpec fileID dataArray ans;

%% SOHR time
filename = fullfile(file_dir, 'output', 'time_SOHR_fraction_all.csv');
delimiter = '';
formatSpec = '%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
fclose(fileID);
SOHR_time = dataArray{:, 1};
clearvars filename delimiter formatSpec fileID dataArray ans;

%% dlsode concentration
filename = fullfile(file_dir, 'output', 'concentration_dlsode_fraction.csv');
delimiter = ',';
formatSpec = '%f%f%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
fclose(fileID);
dlsode_concentration = [dataArray{1:end-1}];
clearvars filename delimiter formatSpec fileID dataArray ans;

%%
labels = {'$O_2$', '$H_2O$', '$H_2$', '$H_2O_2$', '$H$', '$OH$', '$HO_2$', '$O$'};
spe_names = {'O2', 'H2O', 'H2', 'H2O2', 'H', 'OH', 'HO2', 'O'};
% ind = 3;

%% dlsode and figure settings
fig = figure();
pe=plot(dlsode_time, dlsode_concentration(:,ind)-dlsode_concentration(1,ind), 'LineWidth', 2, 'color','r');
hold on;

%% plot SOHR
% load SOHR concentration
filename = fullfile(file_dir, 'output', strcat('concentration_SOHR_fraction_all' ,'.csv'));
delimiter = ',';
formatSpec = '%f%f%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
fclose(fileID);
SOHR_concentration = [dataArray{1:end-1}];
clearvars filename delimiter formatSpec fileID dataArray ans;
% plot
H=plot(SOHR_time, SOHR_concentration(:,ind)-SOHR_concentration(1,ind), '-.', 'LineWidth', 1.5);
hold on;

%% settings
grid on;
%xlim([0, 1.0e-7]);
xlim([dlsode_time(1), dlsode_time(end)]);
set(gca,'GridLineStyle','--');
xlabel('Time/s', 'FontSize', 20);
ylabel('Concentration', 'FontSize', 20);
title(labels(ind), 'FontSize', 20, 'Interpreter','latex');
% use lambda expression
legend_str = 'SOHR';
h_legend = legend([H; pe],legend_str,'EXACT');
set(h_legend, 'FontSize', 14, 'Box', 'off');
set(h_legend, 'Location', 'NorthEast')
if (dlsode_concentration(end, ind)-dlsode_concentration(1,ind))>0
	set(h_legend, 'Location', 'NorthWest')
else
	set(h_legend, 'Location', 'SouthWest')
end

yt = get(gca, 'ytick');
yaxis = get(gca, 'yaxis');
% get the exponent of yticks
expVal = double(yaxis.Exponent);

if expVal ~= 0
    yt_label = strsplit(sprintf('%.1f,',yt/10^expVal), ',');
    set(gca, 'yticklabel', yt_label);  
end

pos = get(gca, 'Position');
x_offset = pos(3)*0.0;
y_offset = - pos(4)*0.015;
if dlsode_concentration(1,ind)>0
    append_str = ['+',num2str(dlsode_concentration(1,ind))];
else
    append_str='';
end

if expVal == 0
    annotation('textbox',[pos(1)+x_offset, pos(2)+pos(4)+y_offset, 0.2, 0.2],...
    'String', append_str,...
    'Interpreter','latex',...
    'VerticalAlignment','bottom',...
    'EdgeColor','none')
else
    annotation('textbox',[pos(1)+x_offset, pos(2)+pos(4)+y_offset, 0.2, 0.2],...
    'String',['$\times10^', '{', num2str(expVal), '}', append_str, '$'],...
    'Interpreter','latex',...
    'VerticalAlignment','bottom',...
    'EdgeColor','none')
end

%% save to file
fname = strjoin(strcat(spe_names(ind), '.png'));
print(fig, fullfile(file_dir, 'output', fname), '-r200', '-dpng');
