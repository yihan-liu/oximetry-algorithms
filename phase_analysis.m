close all
clear;
clc;

%% READ IN
load('cycling_2.mat'); % load data.

TS = 1e-3; % interval for counting

FS = 1 / TS;
START = 5 * FS;

channel_1 = 1;
channel_2 = 3;

s1 = NIR(:, channel_1);
s2 = NIR(:, channel_2);

%% FIRST SMOOTH

hw_s = 0.1; % half-window width in second
hw_sa = hw_s * FS; % half-window width in sample

% initialize the averaged arrays
s1_ave_1 = zeros(length(s1), 1);
s2_ave_1 = zeros(length(s2), 1);

for i =  - hw_sa: hw_sa
    
    s1_ave_1 = s1_ave_1 + circshift(s1, i);
    s2_ave_1 = s2_ave_1 + circshift(s2, i);
    
end

s1_ave_1 = s1_ave_1 ./ (2 * hw_sa + 1);
s2_ave_1 = s2_ave_1 ./ (2 * hw_sa + 1);

% figure
% plot(TIME, s1, TIME, s2);
% hold on

%% PEAK EXTRACTION

[s1_uppper, s1_lower] = envelope(s1_ave_1, 0.8 * FS, 'peak');
[s2_uppper, s2_lower] = envelope(s2_ave_1, 0.8 * FS, 'peak');

p1 = s1_ave_1 - s1_lower;
p2 = s2_ave_1 - s2_lower;

% figure
% plot(TIME, s1, TIME, s1_lower);
% xlim([40 80])
% ylim([480 500])
% hold off

%% SECOND SMOOTH

hw_s = 0.02; % half-window width in second
hw_sa = hw_s * FS; % half-window width in sample

% initialize the averaged arrays
p1_ave = zeros(length(s1), 1);
p2_ave = zeros(length(s2), 1);

for i =  - hw_sa: hw_sa
    
    p1_ave = p1_ave + circshift(p1, i);
    p2_ave = p2_ave + circshift(p2, i);
    
end

p1_ave = p1_ave ./ (2 * hw_sa + 1);
p2_ave = p2_ave ./ (2 * hw_sa + 1);

figure;
plot(TIME, p1, 'g', TIME, p1_ave, 'b',...
    TIME, p2, 'm', TIME, p2_ave, 'r');
hold on

%% FINDPEAKS

[pks1, locs1] = findpeaks(p1_ave, 'MinPeakDistance', 0.4 * FS);
[pks2, locs2] = findpeaks(p2_ave, 'MinPeakDistance', 0.4 * FS);

locs1_s = locs1 ./ FS;
locs2_s = locs2 ./ FS;

plot(locs1_s, pks1, 'ko', locs2_s, pks2, 'kx');
% xlim([44.5, 46.5]);
% ylim([-0.1, 0.1]);

hold on

%% CALCULATE PHASE

locs1_all = [];
locs2_all = [];

len = min(length(locs1_s), length(locs2_s));

for i = 1: len
    
    if min(abs(locs2_s-locs1_s(i))) <= 0.25
        
        locs1_all = [locs1_all; locs1_s(i)];
        
    end
    
    if min(abs(locs1_s-locs2_s(i))) <= 0.25
       
        locs2_all = [locs2_all; locs2_s(i)];
        
    end
    
end

% delete the rightmost element
% locs1_all_s(end) = [];

% these arrays are created only for plots
pks1_all = zeros(length(locs1_all), 1);
pks2_all = zeros(length(locs2_all), 1);

plot(locs1_all, pks1_all, 'ro', locs2_all, pks2_all, 'rx');
% xlim([44.5, 46.5]);
% ylim([-0.1, 0.1]);

hold off

%% VISUALIZE RESULTS

figure

len = min(length(locs1_all), length(locs2_all));
phase_diff = locs2_all(1:len) - locs1_all(1:len);

phase_diff_ave = movmean(phase_diff, 200);

plot(locs2_all(1:len), phase_diff,...
    locs2_all(1:len), abs(phase_diff_ave));
ylim([-0.01, 0.02])

%% EXPORT
% 
% plot_res = [TIME, p1_ave, p2_ave];
% peak_res = [locs1_all_s, locs2_all_s];
% diff_res = [locs2_all_s, phase_diff_ave];

% writematrix(plot_res, 'phase_analysis_plot.csv');
% writematrix(peak_res, 'phase_analysis_peak.csv');
% writematrix(diff_res, 'phase_analysis_diff.csv');
% save('phase_analysis.mat', 'plot_res', 'peak_res', 'diff_res'); 