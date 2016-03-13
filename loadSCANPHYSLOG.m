%
% [OUTPARAMS, WAVEFORMS] = LOADSCANPHYSLOG(SCANPHYSLOG_FILE_NAME)
%
% Loads data from the Philips MRI scanner SCANPHYSLOG file 
% SCANPHYSLOG_FILE_NAME and returns the structures OUTPARAMS and WAVEFORMS
%
% outParams.scanphyslog_file_name % file name of loaded scanphyslog
% outParams.site_name             % site name
% outParams.SRN                   % scanner registration number
% outParams.release               % software release  
% outParams.SWID                  % software identifier
% outParams.date                  % date (yyyy-mm-dd)
% outParams.time                  % time (HH:MM:SS)
% outParams.dockable_table        % flag for dockable table
% outParams.nSamples              % number of samples
% outParams.sample_time_seconds   % sample time spacing
% outParams.time_duration_file_seconds % time duration of the file
% outParams.time_duration_scan_seconds % time duration of the scan
% outParams.scan_start_sample          % sample index of the scan start (mark==0x10)
% outParams.scan_stop_sample           % sample index of the scan stop  (mark==0x20)
%
% waveforms.v1raw        % VKG v1raw waveform
% waveforms.v2raw        % VKG v2raw waveform
% waveforms.v1           % VKG v1 waveform
% waveforms.v2           % VKG v2 waveform
% waveforms.ppu          % peripheral pulse unit (PPU) waveform
% waveforms.resp         % respiratory bellows waverform
% waveforms.gx           % x gradient waveform
% waveforms.gy           % y gradient waveform
% waveforms.gz           % z gradient waveform
% waveforms.mark_hex     % acquisition mark labels hexadecimal number
% waveforms.mark_dec     % acquisition mark labels as decimal number
% waveforms.time_seconds % time in seconds
%
% mark is stored in the file as a hexadecimal bitmask but returned as a
% decimal value. The meaning of the mark values are as follows:
%
% 0x01 = VKG trigger point
% 0x02 = ppu trigger point
% 0x04 = resp trigger point
% 0x08 = measurement marker
% 0x10 = start scan marker
% 0x20 = stop scan marker
%
% Since physiology data is collected from the start of the preparation
% phases, working backward from the stop scan marker by a known scan
% duration is one way to examine the physiology measurements occuring
% during the scan itself. 
%

%
% History:
% 2016.03.10 - welcheb - initial version
% 2016.03.10 - welcheb - version v0.1 for public dissemination
%
function [outParams, waveforms] = loadSCANPHYSLOG(scanphyslog_file_name)

%% initialize outParams
outParams = [];
outParams.scanphyslog_file_name = scanphyslog_file_name;
outParams.site_name = [];
outParams.SRN       = [];
outParams.release   = [];
outParams.SWID      = [];
outParams.date      = [];
outParams.time      = [];
outParams.dockable_table = [];
outParams.nSamples       = [];
outParams.sample_time_seconds         = [];
outParams.time_duration_file_seconds  = [];
outParams.time_duration_scan_seconds  = [];
outParams.scan_start_sample           = [];
outParams.scan_stop_sample            = [];

%% initialize waveforms
waveforms = [];

%% load waveform data
fid = fopen(scanphyslog_file_name);
data = textscan(fid,'%n %n %n %n %n %n %n %n %n %n ','CommentStyle','#');
fclose(fid);
waveforms.v1raw = data{1};
waveforms.v2raw = data{2};
waveforms.v1    = data{3};
waveforms.v2    = data{4};
waveforms.ppu   = data{5};
waveforms.resp  = data{6};
waveforms.gx    = data{7};
waveforms.gy    = data{8};
waveforms.gz    = data{9};
waveforms.mark_hex = data{10};
waveforms.mark_dec = hex2dec( num2str(waveforms.mark_hex) );

%% outParams derived from waveforms
outParams.sample_time_seconds = 0.002; % two milliseconds
outParams.nSamples = numel(waveforms.v1raw);
outParams.scan_start_sample = find(waveforms.mark_dec>=hex2dec('10') & waveforms.mark_dec<hex2dec('20') ); 
outParams.scan_stop_sample  = find(waveforms.mark_dec>=hex2dec('20') ); 
outParams.time_duration_file_seconds = outParams.nSamples * outParams.sample_time_seconds;
outParams.time_duration_scan_seconds = (outParams.scan_stop_sample - outParams.scan_start_sample) * outParams.sample_time_seconds;

%% create time waveform
waveforms.time_seconds = [0:(outParams.nSamples-1)]' * outParams.sample_time_seconds;

%% other outParams information
fid = fopen(scanphyslog_file_name);

%% line 1
str = fgetl(fid);
pat = '##\s+(?<site_name>\w+)\s+\(SRN\s+=\s+(?<SRN>\w+)\),\s+Release\s+(?<release>\w+)\s+\(SWID\s+(?<SWID>\w+)\)';
toks = regexp(str, pat, 'names');
outParams.site_name = toks.site_name;
outParams.SRN       = toks.SRN;
outParams.release   = toks.release;
outParams.SWID      = toks.SWID;

%% line 2
str = fgetl(fid);
pat = '##\s+(?<weekday>\w+)\s+(?<dd>\d+)\-(?<mm>\d+)\-(?<yyyy>\d+)\s+(?<HH>\d+)\:(?<MM>\d+)\:(?<SS>\d+)';
toks = regexp(str, pat, 'names');

outParams.date = sprintf('%s-%s-%s', toks.yyyy, toks.mm, toks.dd );
outParams.time = sprintf('%s:%s:%s', toks.HH, toks.MM, toks.SS );

%% line 3
str = fgetl(fid); % not used

%% line 4
str = fgetl(fid);
pat = '##\s+Dockable\s+table\s+=\s+(?<dockable_table>\w+)';
toks = regexp(str, pat, 'names');
if strcmp(toks.dockable_table,'FALSE'),
    outParams.dockable_table = false;
else
    outParams.dockable_table = true;
end

%% close file
fclose(fid);