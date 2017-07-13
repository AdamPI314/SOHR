%% Import data from text file.
% Script for importing data from the following text file:
%
%    /home/invictus/Documents/spring_2015_boulder/DATA_WORKDIR/TDDM_var_alpha_Mar_30_10_58_23_2015_PP_merge/output/H2O_pathway_correlation_SI/H2O_1st_covariance_SI_topN_pairs_matrix.csv
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2015/04/06 10:03:59

%% Initialize variables.
filename = '/home/invictus/Documents/spring_2015_boulder/DATA_WORKDIR/TDDM_var_alpha_Mar_30_10_58_23_2015_PP_merge/output/H2O_pathway_correlation_SI/H2O_1st_covariance_SI_topN_pairs_matrix.csv';
delimiter = ',';

%% Format string for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
%   column9: double (%f)
%	column10: double (%f)
%   column11: double (%f)
%	column12: double (%f)
%   column13: double (%f)
%	column14: double (%f)
%   column15: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
H2O1stcovarianceSItopNpairsmatrix = [dataArray{1:end-1}];
%% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans;

%% plot 2D bar
fig= figure();
b=bar3(transpose(H2O1stcovarianceSItopNpairsmatrix));
colorbar;

% set color of each surface
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end

% size of matrix
[lx, ly]= size(transpose(H2O1stcovarianceSItopNpairsmatrix));
axis([0.5 ly+.5 0.5 lx+.5 -0.1 0.05]);
xlabel('reaction index');
ylabel('path pair index');
zlabel('$$ 1^{st}$$ order sensitivity index', 'Interpreter', 'Latex');
title('off-diagonal terms');

% sum_SI= sum(sum(H2O1stcovarianceSItopNpairsmatrix) );
% % str = strcat('SUM $$1^{st}$$ order SI: ', num2str(sum_SI));
% str =['SUM $$1^{st}$$ order SI: ', num2str(sum_SI)];
% text(6, 8, 0.10 , str,'Interpreter','latex');

% Change the x and y axis tick labels
xtickpos = 1:1:ly;
ytickpos = 1:1:lx;
% xticklabel= 0:1:ly-1;
xticklabel= [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19];
yticklabel= 0:1:lx-1;
set(gca, 'xtick', xtickpos);
set(gca, 'xticklabel', xticklabel);
set(gca, 'ytick', ytickpos);
set(gca, 'yticklabel', yticklabel);

print(fig, '/home/invictus/Documents/spring_2015_boulder/DATA_WORKDIR/TDDM_var_alpha_Mar_30_10_58_23_2015_PP_merge/output/good_figures/H2O_SI_1st_PATHWAY_correlation_topN_pairs.png', '-r200', '-dpng')