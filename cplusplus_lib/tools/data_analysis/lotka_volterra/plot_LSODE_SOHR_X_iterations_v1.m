% function plot_LSODE_SOHR_X_iterations_v1()

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

%% load sohr concentration
iteration_ind = [1, 3, 5, 7, 9, 11];

concentration_sohr = cell(length(iteration_ind), 1);
counter_iter = 1;
for i_ind=iteration_ind
    % disp(ind)
    filename = fullfile(file_dir, 'output', strcat('concentration_SOHR_M_all_', num2str(i_ind) ,'.csv'));
    delimiter = ',';
    formatSpec = '%f%f%f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
    fclose(fileID);
    concentration_sohr{counter_iter, 1} = [dataArray{1:end-1}];
    counter_iter = counter_iter + 1;
    clearvars filename delimiter formatSpec fileID dataArray ans;
end

%%
% colors = {'r','b','k','c','m','y','g','w'};
% colors = hsv(length(iteration_ind) + 1);
colors = jet(length(iteration_ind) + 1);
marker = {'hexagram', 'o', '*', '.', 'x', 'square', 'diamond', '^', 'v', '>', '<', 'pentagram', '+'};


%% dlsode and figure settings
fig = figure();
ratio = 1.0/3.0;
species_ind=2;

Nexact = 10;
Nsohr=1;

% Handler array
H = gobjects(length(iteration_ind) + 1);
color_counter = 1;
marker_counter = 1;
hanlder_counter = 1;

%% plot sohr
for i_ind=1:length(iteration_ind)    
    H(hanlder_counter) = plot(sohr_time(1:Nsohr:floor(end*ratio)), concentration_sohr{i_ind,1}(1:Nsohr:floor(end*ratio), species_ind), '-', 'LineWidth', 1.5, 'color', colors(color_counter, :));
    color_counter = color_counter + 1;
    hanlder_counter = hanlder_counter + 1;
    hold on;
end

%% plot lsode
H(hanlder_counter) = scatter(dlsode_time(1:Nexact:floor(end*ratio)), dlsode_concentration(1:Nexact:floor(end*ratio),species_ind), 'marker', marker{1, marker_counter}, 'MarkerFaceColor', colors(color_counter, :), 'MarkerEdgeColor', colors(color_counter, :));
hanlder_counter = hanlder_counter + 1;
hold on;

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

% xlim([dlsode_time(1), dlsode_time(floor(end*ratio))]);

xlabel('Time', 'FontSize', 20);
ylabel('Concentration', 'FontSize', 20);

% use lambda expression
legend_str = cell(1, 1+length(iteration_ind));
legend_str{1, end} = 'EXACT';
for iter_ind = 1:length(iteration_ind)
    legend_str{1, iter_ind} = ['n = ', num2str(iteration_ind(iter_ind))];
end 

h_legend = legend(H(:,1) ,legend_str);
set(h_legend, 'FontSize', 16, 'Box', 'off');
set(h_legend, 'Location', 'NorthEast');

h_pos=get(h_legend,'position');
set(h_legend, 'position', [h_pos(1)*0.7, h_pos(2)*0.8, h_pos(3), h_pos(4)]);

%% save to file
fname = strcat('LSODE_SOHR_X_iterations', '.png');
print(fig, fullfile(file_dir, 'output', fname), '-r200', '-dpng');
