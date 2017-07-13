%% Import data from text file.
% Script for importing data from the following text file:
%
%    /home/invictus/Documents/spring_2015_boulder/DATA_WORKDIR/TDDM_var_alpha_Mar_30_10_58_23_2015_PP_merge/output/H2O_pathway_correlation_SI/H2O_1st_covariance_SI_matrix.csv
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2015/04/02 15:51:30

%% Initialize variables.
filename = '/home/invictus/Documents/spring_2015_boulder/DATA_WORKDIR/TDDM_var_alpha_Mar_30_10_58_23_2015_PP_merge/output/H2O_pathway_correlation_SI/H2O_1st_covariance_SI_matrix.csv';
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
%	column16: double (%f)
%   column17: double (%f)
%	column18: double (%f)
%   column19: double (%f)
%	column20: double (%f)
%   column21: double (%f)
%	column22: double (%f)
%   column23: double (%f)
%	column24: double (%f)
%   column25: double (%f)
%	column26: double (%f)
%   column27: double (%f)
%	column28: double (%f)
%   column29: double (%f)
%	column30: double (%f)
%   column31: double (%f)
%	column32: double (%f)
%   column33: double (%f)
%	column34: double (%f)
%   column35: double (%f)
%	column36: double (%f)
%   column37: double (%f)
%	column38: double (%f)
%   column39: double (%f)
%	column40: double (%f)
%   column41: double (%f)
%	column42: double (%f)
%   column43: double (%f)
%	column44: double (%f)
%   column45: double (%f)
%	column46: double (%f)
%   column47: double (%f)
%	column48: double (%f)
%   column49: double (%f)
%	column50: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

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
H2O1stcovarianceSImatrix = [dataArray{1:end-1}];
%% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans;

%% plot 2D bar
fig= figure();
b=bar3(H2O1stcovarianceSImatrix(1:10,1:10));
colorbar;

% set color of each surface
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end

axis([0.5 10.5 0.5 10.5 -0.1 0.05]);
xlabel('path index');
ylabel('path index');
zlabel('$$ 1^{st}$$ order sensitivity index', 'Interpreter', 'Latex');
title('off-diagonal terms');

sum_SI= sum(sum(H2O1stcovarianceSImatrix(1:10,1:10)) );
% str = strcat('SUM $$1^{st}$$ order SI: ', num2str(sum_SI));
str =['SUM $$1^{st}$$ order SI: ', num2str(sum_SI)];
text(6, 8, 0.10 , str,'Interpreter','latex');

% Change the x and y axis tick labels
xtickpos = 1:1:10;
ytickpos = 1:1:10;
xticklabel= 0:1:9;
yticklabel= 0:1:9;
set(gca, 'xtick', xtickpos);
set(gca, 'xticklabel', xticklabel);
set(gca, 'ytick', ytickpos);
set(gca, 'yticklabel', yticklabel);

print(fig, '/home/invictus/Documents/spring_2015_boulder/DATA_WORKDIR/TDDM_var_alpha_Mar_30_10_58_23_2015_PP_merge/output/good_figures/H2O_SI_1st_PATHWAY_correlation_all.png', '-r200', '-dpng')