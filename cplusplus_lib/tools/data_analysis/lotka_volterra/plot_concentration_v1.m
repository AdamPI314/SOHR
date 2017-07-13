% function plot_concentration_v2(ind)
% ind =8;

%% Current file directory
file_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', '..');

%% dlsode time
filename = fullfile(file_dir, 'output', 'time_dlsode_M.csv');
delimiter = '';
formatSpec = '%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
fclose(fileID);
dlsode_time = dataArray{:, 1};
clearvars filename delimiter formatSpec fileID dataArray ans;

%% sohr time
filename = fullfile(file_dir, 'output', 'time_SOHR_M_all.csv');
delimiter = '';
formatSpec = '%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
fclose(fileID);
sohr_time = dataArray{:, 1};
clearvars filename delimiter formatSpec fileID dataArray ans;

%% dlsode concentration
filename = fullfile(file_dir, 'output', 'concentration_dlsode_M.csv');
delimiter = ',';
formatSpec = '%f%f%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
fclose(fileID);
dlsode_concentration = [dataArray{1:end-1}];
clearvars filename delimiter formatSpec fileID dataArray ans;

%%
labels = {'$X$', '$Y$'};
spe_names = {'X', 'Y'};
colors = {'r','b','k','c','m','y','g','w'};
marker = {'hexagram', 'o', '*', '.', 'x', 'square', 'diamond', '^', 'v', '>', '<', 'pentagram', '+'};
%% dlsode and figure settings
fig = figure();
counter = 1;
Nexact = 5;
% Handler array
H = gobjects(4);
h_counter = 1;

for ind=2:3
    H(h_counter) = scatter(dlsode_time(1:Nexact:end), dlsode_concentration(1:Nexact:end,ind), 'marker', marker{1, counter}, 'MarkerFaceColor', colors{1,counter}, 'MarkerEdgeColor', colors{1,counter});
    counter = counter + 1;
    h_counter = h_counter + 1;
    hold on;
end

%% plot sohr
% load sohr concentration
filename = fullfile(file_dir, 'output', strcat('concentration_sohr_M_all' ,'.csv'));
delimiter = ',';
formatSpec = '%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
fclose(fileID);
sohr_concentration = [dataArray{1:end-1}];
clearvars filename delimiter formatSpec fileID dataArray ans;

% plot
counter = 1;
Nsohr=100;
for ind=2:3
    H(h_counter) = plot(sohr_time(1:Nsohr:end), sohr_concentration(1:Nsohr:end, ind), '-.', 'LineWidth', 1.5, 'color', colors{1,counter});
    counter = counter + 1;
    h_counter = h_counter + 1;
    hold on;
end

%% settings
grid on;
%xlim([0, 1.0e-7]);

xticks = get(gca, 'xtick');
set(gca, 'xtick', xticks);
xticklabel = get(gca, 'XTickLabel');
set(gca,'XTickLabel', xticklabel,'FontName','Times','fontsize',14);
yticks = get(gca, 'ytick');
set(gca, 'ytick', yticks);
yticklabel = get(gca, 'YTickLabel');
set(gca,'YTickLabel', yticklabel,'FontName','Times','fontsize',14);

xlim([dlsode_time(1), dlsode_time(end)]);
xlabel('Time', 'FontSize', 20);
ylabel('Concentration', 'FontSize', 20);
title('Concentrations of intermediates X and Y', 'FontSize', 20, 'Interpreter','latex');

% use lambda expression
legend_str = {'EXACT X', 'EXACT Y', 'SOHR X', 'SOHR Y'};
h_legend = legend([H(1), H(2), H(3), H(4)],legend_str);
set(h_legend, 'FontSize', 12, 'Box', 'off');
set(h_legend, 'Location', 'NorthWest');

h_pos=get(h_legend,'position');
set(h_legend, 'position', [h_pos(1)*1.35, h_pos(2), h_pos(3), h_pos(4)]);

%% save to file
fname = strcat('LSODE_SOHR_X_Y', '.png');
print(fig, fullfile(file_dir, 'output', fname), '-r200', '-dpng');
