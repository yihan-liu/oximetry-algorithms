close all;
clearvars -except data file_name;
clc;

%% READ IN

if ~exist('data', 'var')
    
    file_name = 'datasets/cycling_4';
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
TS = 1.2e-3; % (s/sample) interval between two neighboring selected points
channel = 2; % choose a channel to calculate

RED_offset = 0.7e-3; % first point of RED signal
NIR_offset = 0.1e-3; % first point of NIR signal

% plot(data(:, 1), data(:, 2)) % for test

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
    FrequencyLimits=[0.5, 3]);
[wt, f] = cwt(RED(:, channel), FilterBank=fb);

figure
cwt(RED(:, channel), FilterBank=fb);

colorbar off

wt_ave = movmean(abs(wt), 5*FS, 2);

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

%% HEART RATE EXTRACTION

wtt = wt_ave(f>0.1&f<2, :);
cut = length(f(f>=2));

[~, freq_index] = max(wtt, [], 1);
max_freq = f(freq_index+cut);

figure
plot(TIME, max_freq);
% writematrix([TIME, RED(:, channel), NIR(:, channel)], 'raw_data_hypoxia6_channel3.csv');
% writematrix([TIME, max_freq], 'heart_rate_res_hypoxia6_channel3.csv');