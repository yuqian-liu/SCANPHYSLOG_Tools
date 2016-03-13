%% Run example SCANPHYSLOG processing

%% start with a clean slate
close all; clear all; clc;

%% parent folder name
parent_folder_name = './data_input';

%% processing parameters
processParams.heart_rate_time_window_seconds           = 30;
processParams.heart_rate_time_window_overlap_fraction  = 0.5;
processParams.breath_rate_time_window_seconds          = 60;
processParams.breath_rate_time_window_overlap_fraction = 0.5;
processParams.confidence_level_interval_probability    = 0.95;
processParams.heart_rate_max  = 200;
processParams.heart_rate_min  = 30;
processParams.breath_rate_max = 30;
processParams.breath_rate_min = 8;
processParams.match_pattern       = '*SCANPHYSLOG*';
processParams.output_csv_filename = 'SCANPHYSLOG.csv';
processParams.output_folder = './data_output';
processParams.verbose = true;

%% results
outParams = process_SCANPHYSLOG_parent_folder(parent_folder_name, processParams);

%% load CSV file as a MATLAB Table
T = readtable( sprintf('%s/%s', processParams.output_folder, processParams.output_csv_filename) );

%% display full Table
T

%% display abbreviated summary
summary(T(:,{'heart_rate_mean','breath_rate_mean'}));

%% create plot
figure('Position',[50, 50, 1000, 500]);

subplot(1,2,1);
errorbar(T{:,{'idx'}}, T{:,{'heart_rate_mean'}},  T{:,{'heart_rate_CI_HW'}});
xlabel('File Index');
ylabel('Beats per Minute');
title('Heart Rate');
grid on;

subplot(1,2,2);
errorbar(T{:,{'idx'}}, T{:,{'breath_rate_mean'}},  T{:,{'breath_rate_CI_HW'}});
xlabel('File Index');
ylabel('Breaths per Minute');
title('Respiratory Rate');
grid on;

%% save PNG file
png_filename = sprintf('%s/SCANPHYSLOG.png', processParams.output_folder);
imwrite(frame2im(getframe(gcf)), png_filename);