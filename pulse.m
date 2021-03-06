close all;
clearvars -except data file_name;
clc;

%% READ IN

if ~exist('data', 'var')
    
    file_name = 'datasets/wrist_hypoxia_4';
    format = '.txt';
    fid = fopen(strcat(file_name, format));
    data = cell2mat(textscan(fid, '%f %f %f %f %f',...
        Delimiter='\t',...
        HeaderLines=6));
    fclose(fid);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters that need to be specified every time %

FS_raw = 10000; % (Hz) sampling frequency of raw data
TS = 1e-3; % (s/sample) interval between two neighboring selected points
channel = 4; % choose a channel to calculate

RED_offset = 0.6e-3; % first point of RED signal
NIR_offset = 0.1e-3; % first point of NIR signal

% figure; plot(data(:, 1), data(:, 4)) % for test

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TS_raw = 1/FS_raw; % (s/sample)
TIME_raw = data(:, 1)-data(1, 1); % stores the raw time data
LEN = length(data); % length of data
FS = 1/TS; % (sample/s) new sampling frequency

TIME = transpose(0:TS:TIME_raw(end));
RED = data(ceil(RED_offset*FS_raw):ceil(TS*FS_raw):LEN, 2:5);
NIR = data(ceil(NIR_offset*FS_raw):ceil(TS*FS_raw):LEN, 2:5);

% RED(end+1, :) = 0;

% save(strcat(file_name, '.mat'), 'TIME', 'RED', 'NIR'); 
% save as .mat file

%% WAVELET TRANSFORM

fb = cwtfilterbank(...
    SignalLength=length(TIME),...
    VoicesPerOctave=48,...
    SamplingFrequency=FS,...
    FrequencyLimits=[0.5, 3],...
    TimeBandwidth=120);
[wt, f] = cwt(RED(:, channel), FilterBank=fb);

figure
cwt(RED(:, channel), FilterBank=fb);
caxis([0 0.6]);
title("");
xlabel([]);
colorbar off

wt_ave = movmean(abs(wt), 1*FS, 2);

%% PLOT 
 
% freq_pulse = 68; % position of the pulse frequency
% freq_resp = 100; % position of the respiratory frequency
% 
% % figure;
% % cwt(NIR_1, FS);
% %  
% % figure;
% % cwt(NIR_2, FS);
% 
% figure(1);
% plot(TIME, wt_ave(freq_pulse, :), 'r');
% title(strcat('Pulse magnitude', ' @ ', string(f2(freq_pulse)), 'Hz'));
% ylim([0, 0.1]);
% 
% figure(2);
% plot(TIME, wt_ave(freq_resp, :), 'r');
% title(strcat('Respiratory magnitude', ' @ ', string(f2(freq_resp)), 'Hz'));
% ylim([0, 0.5]);
% 
% pulse_res = [TIME, wt_ave(freq_pulse, :)', wt_ave(freq_resp, :)'];
% % 
% % % writematrix(res, 'pulse & respiratory.csv');
% % save('pulse_extraction.mat', 'pulse_res');

%% DATA EXTRACTION

upper_p = 1.2;
lower_p = 0.8;
upper_r = 0.7;
lower_r = 0.5;

wtt_p = wt_ave(f>lower_p&f<upper_p, :);
cut_p = length(f(f>=upper_p));

wtt_r = wt_ave(f>lower_r&f<upper_r, :);
cut_r = length(f(f>=upper_r));

% HR and RR

[~, hr_freq_index] = max(wtt_p, [], 1);
hr_freq_index = hr_freq_index';
hr_max_freq = f(hr_freq_index+cut_p);
HR = 60*hr_max_freq;

[~, rr_freq_index] = max(wtt_r, [], 1);
rr_freq_index = rr_freq_index';
rr_max_freq = f(rr_freq_index+cut_r);
RR = 60*rr_max_freq;

% PI and RI

% steady
% PI = wtt(find(f>1.1&f<1.11)-cut, :)'; %#ok<FNDSB>
% RI = wtt(find(f>0.5&f<0.505)-cut, :)'; %#ok<FNDSB>

% moving
PI = zeros(length(TIME), 1);
for i = 1: length(TIME)
    PI(i) = wtt_p(hr_freq_index(i), i);
end

RI = zeros(length(TIME), 1);
for i = 1: length(TIME)
    RI(i) = wtt_r(rr_freq_index(i), i);
end

figure
plot(TIME, HR);

% surf(wtt, EdgeColor='none');
% writematrix([TIME, RED(:, channel), NIR(:, channel)], 'raw_data_hypoxia6_channel3.csv');
writematrix([TIME, HR, RR], 'HR_RR_res_hyp4_ch4.csv');