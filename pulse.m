% close all;
clearvars -except data file_name;
clc;

%% READ IN

if ~exist('data', 'var')
    
    file_name = 'datasets/wrist_hypoxia_3';
    format = '.txt';
    fid = fopen(strcat(file_name, format));
    data = cell2mat(textscan(fid, '%f %f %f %f %f',...
        Delimiter='\t',...
        HeaderLines=6));
    fclose(fid);
    
end

FS_raw = 10000; % s
TS_raw = 1/FS_raw; % samples/s
TIME_raw = data(:, 1)-data(1, 1);

channel = 4;

RED_offset = 0.6e-3;
NIR_offset = 0.1e-3;
LEN = length(data);

TS = 1e-3; % s
FS = 1/TS; % samples/s

TIME = transpose(0:TS:TIME_raw(end));
RED = data(ceil(RED_offset*FS_raw):ceil(TS*FS_raw):LEN, 2:5);
NIR = data(ceil(NIR_offset*FS_raw):ceil(TS*FS_raw):LEN, 2:5);

plot(TIME, NIR);
hold on
plot(TIME, RED)
hold off
% plot(TIME_raw, data(:, 2));
% save(strcat(file_name, '.mat'), 'TIME', 'RED', 'NIR');

%% WAVELET TRANSFORM

fb = cwtfilterbank(...
    SignalLength=length(TIME),...
    VoicesPerOctave=48,...
    SamplingFrequency=FS,...
    FrequencyLimits=[0.5, 2]);
[wt, f] = cwt(NIR(:, channel), FilterBank=fb);
% cwt(NIR(:, channel), FS)

wt_ave = movmean(abs(wt), 5 * FS, 2);

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
% 
% % writematrix(res, 'pulse & respiratory.csv');
% save('pulse_extraction.mat', 'pulse_res');

%% HEART RATE EXTRACTION

wtt = wt_ave(f>0.1&f<1.5, :);
cut = length(f(f>=2));

[~, freq_index] = max(wtt, [], 1);
max_freq = f(freq_index+cut);

plot(TIME, max_freq);
% writematrix([TIME, max_freq], 'heart_rate_res_hyp4_channel2.csv');
