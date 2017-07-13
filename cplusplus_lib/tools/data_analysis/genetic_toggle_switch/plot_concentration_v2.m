function plot_concentration_v2()

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
colors = {'r','b','k','c','m','y','g','w'};
Nlines = 4;
% colors = jet(Nlines/2);

marker = {'hexagram', 'o', '*', '.', 'x', 'square', 'diamond', '^', 'v', '>', '<', 'pentagram', '+'};
%% dlsode and figure settings
fig = figure();

counter = 1;
Nexact = 1;
Nsohr=25;

% Handler array
H = gobjects(Nlines);
h_counter = 1;

for ind=2:3
%     H(h_counter) = plot(dlsode_time(1:Nexact:end), dlsode_concentration(1:Nexact:end,ind), '-', 'LineWidth', 1.5, 'color', colors(counter, :));
    H(h_counter) = plot(dlsode_time(1:Nexact:end), dlsode_concentration(1:Nexact:end,ind), '-', 'LineWidth', 1.5, 'color', colors{1, counter});
    counter = counter + 1;
    h_counter = h_counter + 1;
    hold on;
end

%% plot sohr
% load sohr concentration
filename = fullfile(file_dir, 'output', strcat('concentration_SOHR_M_all' ,'.csv'));
delimiter = ',';
formatSpec = '%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
fclose(fileID);
sohr_concentration = [dataArray{1:end-1}];
clearvars filename delimiter formatSpec fileID dataArray ans;

% plot
counter = 1;

for ind=2:3    
%     H(h_counter) = scatter(sohr_time(1:Nsohr:end), sohr_concentration(1:Nsohr:end, ind), 'marker', marker{1, counter}, 'MarkerFaceColor', colors(counter, :), 'MarkerEdgeColor', colors(counter, :));
    H(h_counter) = scatter(sohr_time(1:Nsohr:end), sohr_concentration(1:Nsohr:end, ind), 'marker', marker{1, counter}, 'MarkerFaceColor', colors{1, counter}, 'MarkerEdgeColor', colors{1, counter});
    counter = counter + 1;
    h_counter = h_counter + 1;
    hold on;
end

%% settings
grid on;

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
% title('Concentrations of intermediate X and Y', 'FontSize', 20, 'Interpreter','latex');

% use lambda expression
legend_str = {'EXACT X', 'EXACT Y', 'SOHR X', 'SOHR Y'};
h_legend = legend(H(:,1),legend_str);
set(h_legend, 'FontSize', 12, 'Box', 'off');
set(h_legend, 'Location', 'NorthWest');

h_pos=get(h_legend,'position');
set(h_legend, 'position', [h_pos(1)*1.25, h_pos(2), h_pos(3), h_pos(4)]);

%% save to file
fname = strcat('LSODE_SOHR_X_Y', '.png');
print(fig, fullfile(file_dir, 'output', fname), '-r200', '-dpng');
