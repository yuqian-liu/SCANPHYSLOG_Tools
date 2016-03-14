%
% OUTPARAMS = PROCESS_SCANPHYSLOG_PARENT_FOLDER(PARENT_FOLDER_NAME)
%
% Processes all SCANPHYSLOG files in the folder PARENT_FOLDER_NAME and
% returns the array of structures OUTPARAMS containing
%
% outParams(:).date               % date (yyyy-mm-dd)
% outParams(:).time               % time (HH:MM:SS)    
% outParams(:).N_heart_rate       % number of heart rate windows
% outParams(:).heart_rate_mean    % heart rate mean
% outParams(:).heart_rate_sigma   % heart rate standard deviation
% outParams(:).heart_rate_CI_HW % heart rate confidence interval half-width
% outParams(:).heart_rate_CI_LB % heart rate confidence interval lower bound
% outParams(:).heart_rate_CI_UB % heart rate confidence interval upper bound
% outParams(:).N_breath_rate       % number of breath rate windows
% outParams(:).breath_rate_mean    % breath rate mean
% outParams(:).breath_rate_sigma   % breath rate standard deviation
% outParams(:).breath_rate_CI_HW % breath rate confidence interval half-width
% outParams(:).breath_rate_CI_LB % breath rate confidence interval lower bound
% outParams(:).breath_rate_CI_UB % breath rate confidence interval upper bound
%
% OUTPARAMS = PROCESS_SCANPHYSLOG_PARENT_FOLDER(PARENT_FOLDER_NAME, PROCESSPARAMS)
%
% A comma separated variable file names 'SCANPHYSLOG.csv' is also written
% to disk.
%
% Uses the structure PROCESSPARAMS to control the processing of the
% SCANPHYSLOG files. Options and default values are:
%
% processParams.heart_rate_time_window_seconds           = 30;  % heart rate time window duration in seconds
% processParams.heart_rate_time_window_overlap_fraction  = 0.5; % heart rate window overlap fraction
% processParams.breath_rate_time_window_seconds          = 60;  % breath rate time window duration in seconds
% processParams.breath_rate_time_window_overlap_fraction = 0.5; % breath rate window overlap fraction
% processParams.confidence_level_interval_probability    = 0.95; % level for reported confidence intervals
% processParams.heart_rate_max  = 200;                           % maxmimum detectable heart rate (beats per minute)
% processParams.heart_rate_min  = 30;                            % minimum detectable heart rate (beats per minute)
% processParams.breath_rate_max = 30;                            % maxmimum detectable breath rate (breaths per minute)
% processParams.breath_rate_min = 8;                             % minimum detectable breath rate (breaths per minute)
% processParams.match_pattern       = '*SCANPHYSLOG*';           % pattern to detect SCANPHYSLOG files to process
% processParams.output_csv_filename = 'SCANPHYSLOG.csv';         % name for output CSV file
% processParams.output_folder = PARENT_FOLDER_NAME;              % name of parent folder containing SCANPHYSLOG files
% processParams.verbose = true;                                  % print progress messages
%
% Developed with MATLAB R2014b
%
% Supported by a grant from the National Institutes of Health (NIH), National
% Institute of Diabetes and Digestive and Kidney Diseases (NIDDK) - R01DK105371
%

%
% History:
% 2016.03.01 - welcheb - initial version
% 2016.03.10 - welcheb - version v0.1 for public dissemination
%
function [outParams] = process_SCANPHYSLOG_parent_folder(parent_folder_name, processParams)

%% handle processParams
if nargin<2,
    processParams = [];
end
if ~isfield(processParams,'heart_rate_time_window_seconds'),
    processParams.heart_rate_time_window_seconds = 30;
end
if ~isfield(processParams,'heart_rate_time_window_overlap_fraction'),
    processParams.heart_rate_time_window_overlap_fraction = 0.5;
end
if ~isfield(processParams,'breath_rate_time_window_seconds'),
    processParams.breath_rate_time_window_seconds = 30;
end
if ~isfield(processParams,'breath_rate_time_window_overlap_fraction'),
    processParams.breath_rate_time_window_overlap_fraction = 0.5;
end
if ~isfield(processParams,'confidence_level_probability'),
    processParams.confidence_level_probability = 0.95;
end
% Normal resting heart rate for healthy human adult is 60 to 100 beats per minute
if ~isfield(processParams,'heart_rate_max'),
    processParams.heart_rate_max = 200;
end
if ~isfield(processParams,'heart_rate_min'),
    processParams.heart_rate_min = 30;
end
% Normal resting breath rate for healthy human adult is 12 t 25 breaths per minute
if ~isfield(processParams,'breath_rate_max'),
    processParams.breath_rate_max = 30;
end
if ~isfield(processParams,'breath_rate_min'),
    processParams.breath_rate_min = 8;
end
if ~isfield(processParams,'match_pattern'),
    processParams.match_pattern = '*SCANPHYSLOG*';
end
if ~isfield(processParams,'output_csv_filename'),
    processParams.output_csv_filename = 'SCANPHYSLOG.csv';
end
if ~isfield(processParams,'output_folder'),
    processParams.output_folder = parent_folder_name;
end
if ~isfield(processParams,'verbose'),
    processParams.verbose = true;
end

%% initialize outParams
outParams = [];
outParams(1).date = [];
outParams(1).time = [];
outParams(1).N_heart_rate= [];
outParams(1).heart_rate_mean = [];
outParams(1).heart_rate_sigma = [];
outParams(1).heart_rate_CI_HW = [];
outParams(1).heart_rate_CI_LB = [];
outParams(1).heart_rate_CI_UB = [];
outParams(1).N_breath_rate = [];
outParams(1).breath_rate_mean = [];
outParams(1).breath_rate_sigma = [];
outParams(1).breath_rate_CI_HW = [];
outParams(1).breath_rate_CI_LB = [];
outParams(1).breath_rate_CI_UB = [];

%% get list of scanphyslog files
d = dir( sprintf('%s/%s', parent_folder_name, processParams.match_pattern) );

%% setup power spectrum calculation
file_sample_period_seconds = 2e-3; % 2 ms
fs = 1/file_sample_period_seconds;

%% heart rate
heart_rate_segmentLength = round( processParams.heart_rate_time_window_seconds / file_sample_period_seconds );
heart_rate_nOverlap = round(heart_rate_segmentLength * processParams.heart_rate_time_window_overlap_fraction);

%% breath rate
breath_rate_segmentLength = round( processParams.breath_rate_time_window_seconds / file_sample_period_seconds );
breath_rate_nOverlap = round(breath_rate_segmentLength * processParams.breath_rate_time_window_overlap_fraction);

%% initialize outParams
outParams = [];

%% loop over files
for idx = 1:numel(d),
    
    %% load SCANPHYSLOG file
    [out, waveforms] = loadSCANPHYSLOG( sprintf('%s/%s', parent_folder_name, d(idx).name) );
    PPU = waveforms.ppu;
    resp = waveforms.resp;
    
    outParams(idx).date = out.date;
    outParams(idx).time = out.time;
    
    %% butterworth bandpass filter PPU 
	[PPU_b, PPU_a] = butter(3, [processParams.heart_rate_min/60 processParams.heart_rate_max/60]/fs);
    PPU_filt = filter(PPU_b, PPU_a, PPU - mean(PPU));
    
    %% butterworth bandpass filter resp 
	[resp_b, resp_a] = butter(3, [processParams.breath_rate_min/60 100*processParams.breath_rate_max/60]/fs);
    resp_filt = filter(resp_b, resp_a, resp - mean(resp) );
    
    %% loop for overlapping PPU windows
    window_start = 1;
    window_stop  = window_start + heart_rate_segmentLength - 1;
	delta_f_PPU = 1/( heart_rate_segmentLength * out.sample_time_seconds );
    if mod(heart_rate_segmentLength,2)==0;
        f_PPU_v = [-heart_rate_segmentLength/2:(+heart_rate_segmentLength/2-1)] * delta_f_PPU;
    else
        f_PPU_v = [-(heart_rate_segmentLength-1)/2:+(heart_rate_segmentLength-1)/2] * delta_f_PPU;
    end
    f_PPU_v_idx_examine = find( f_PPU_v>=(processParams.heart_rate_min/60) & f_PPU_v<=(processParams.heart_rate_max/60) );
    f_PPU_examine = f_PPU_v(f_PPU_v_idx_examine);
    PPU_heart_beats_per_minute_v = [];
    k = 1;
    while window_stop < out.nSamples,
        fft_PPU = abs(fftshift(fft(PPU_filt([window_start:window_stop]))));
        fft_PPU_examine = fft_PPU(f_PPU_v_idx_examine);
        peak_PPU_examine_idx = find(abs(fft_PPU_examine)==max(abs(fft_PPU_examine)));
        peak_PPU_examine_idx = peak_PPU_examine_idx(1);
        PPU_heart_beats_per_minute_v(k) = 60 * f_PPU_examine(peak_PPU_examine_idx);
        window_start = window_start + heart_rate_segmentLength;
        window_stop = window_stop + heart_rate_segmentLength;
        k = k + 1;
    end
    
    %% loop for overlapping resp windows
    window_start = 1;
    window_stop  = window_start + breath_rate_segmentLength - 1;
	delta_f_resp = 1/( breath_rate_segmentLength * out.sample_time_seconds );
    if mod(breath_rate_segmentLength,2)==0;
        f_resp_v = [-breath_rate_segmentLength/2:(+breath_rate_segmentLength/2-1)] * delta_f_resp;
    else
        f_resp_v = [-(breath_rate_segmentLength-1)/2:+(breath_rate_segmentLength-1)/2] * delta_f_resp;
    end
    f_resp_v_idx_examine = find( f_resp_v>=(processParams.breath_rate_min/60) & f_resp_v<=(processParams.breath_rate_max/60) );
    f_resp_v_examine = f_resp_v(f_resp_v_idx_examine);
    resp_breaths_per_minute_v = [];
    k = 1;
    while window_stop < out.nSamples,
        fft_resp = abs(fftshift(fft(resp_filt([window_start:window_stop]))));
        fft_resp_examine = fft_resp(f_resp_v_idx_examine);
        peak_resp_examine_idx = find(abs(fft_resp_examine)==max(abs(fft_resp_examine)));
        peak_resp_examine_idx = peak_resp_examine_idx(1);
        resp_breaths_per_minute_v(k) = 60 * f_resp_v_examine(peak_resp_examine_idx);
        window_start = window_start + breath_rate_segmentLength;
        window_stop = window_stop + breath_rate_segmentLength;
        k = k + 1;
    end
        
    %% other entries of outParams
    
    % heart rate
    outParams(idx).N_heart_rate = numel(PPU_heart_beats_per_minute_v);
    outParams(idx).heart_rate_mean = mean(PPU_heart_beats_per_minute_v);
    outParams(idx).heart_rate_sigma = std(PPU_heart_beats_per_minute_v);
    tval_heart_rate = tinv( processParams.confidence_level_interval_probability, outParams(idx).N_heart_rate - 1);
    outParams(idx).heart_rate_CI_HW  = tval_heart_rate * outParams(idx).heart_rate_sigma / sqrt( outParams(idx).N_heart_rate );
    outParams(idx).heart_rate_CI_LB = outParams(idx).heart_rate_mean - outParams(idx).heart_rate_CI_HW;
    outParams(idx).heart_rate_CI_UB = outParams(idx).heart_rate_mean + outParams(idx).heart_rate_CI_HW;
    
	% breath rate
    outParams(idx).N_breath_rate = numel(resp_breaths_per_minute_v);
    outParams(idx).breath_rate_mean = mean(resp_breaths_per_minute_v);
    outParams(idx).breath_rate_sigma = std(resp_breaths_per_minute_v);
    tval_breath_rate = tinv( processParams.confidence_level_interval_probability, outParams(idx).N_breath_rate - 1);
    outParams(idx).breath_rate_CI_HW = tval_breath_rate * outParams(idx).breath_rate_sigma / sqrt( outParams(idx).N_breath_rate );
    outParams(idx).breath_rate_CI_LB = outParams(idx).breath_rate_mean - outParams(idx).breath_rate_CI_HW;
    outParams(idx).breath_rate_CI_UB = outParams(idx).breath_rate_mean + outParams(idx).breath_rate_CI_HW;
    
    %% verbose message
    if processParams.verbose,
        disp( sprintf('processed SCANPHYSLOG file %03d of %03d', idx, numel(d)) );
    end
    
end

%% save to CSV file
fid = fopen( sprintf('%s/%s', processParams.output_folder, processParams.output_csv_filename), 'w');
outParams_fieldnames = fieldnames(outParams);
fprintf(fid,'"idx"');
for k = 1:numel(outParams_fieldnames),
    fn = outParams_fieldnames{k};
    fprintf(fid,',"%s"', fn);
end
fprintf(fid,'\n');

for idx = 1:numel(d),
    fprintf(fid, '%d', idx);
    for k = 1:numel(outParams_fieldnames),
        fn = outParams_fieldnames{k};
        val = outParams(idx).(fn);
        if ischar(val),
            fprintf(fid, ',"%s"', val);
        else
            fprintf(fid, ',%.2f', val);
        end
    end
    fprintf(fid,'\n');
end