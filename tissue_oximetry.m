% close all
clear;
clc;

%% READ IN
load('datasets/wrist_hypoxia_7.mat');

TS = 1e-3; % interval for counting

FS = 1/TS;
START = 5*FS;

% select channel
channel = 4;
NIR_raw = NIR(:, channel);
RED_raw = RED(:, channel);

% trim data
if length(RED) > length(NIR)
    
    RED_raw(length(NIR)+1:end, :) = [];
    TIME(length(NIR)+1:end, :) = [];
    
elseif length(RED) < length(NIR)
    
    NIR_raw(length(RED_raw)+1:end, :) = [];
    TIME(length(RED_raw)+1:end, :) = [];
    
end

% plot the raw data
% subplot(4, 2, 1);
% plot(TIME(START+1:end),...
%     NIR(START+1:end,channel),...
%     'b');
% title('Raw Data', fontsize=30);
% ylabel('Voltage(mV)', fontsize=16);
% set(gca, Fontsize=14);
% xlim([0, TIME(end)]);
% % ylim([300, 400]);
% 
% subplot(4, 2, 3);
% plot(TIME(START+1:end),...
%     RED(START+1:end, channel),...
%     'r');
% xlabel('Time(s)', fontsize=16);
% ylabel('Voltage(mV)', fontsize=16);
% set(gca, Fontsize=14);
% xlim([0, TIME(end)]);
% xlim([0, TIME(end)]);
% % ylim([350, 450]);



%% DENOISE (LOCAL AVERAGE)

NIR_ave = movmean(NIR_raw, 5*FS, 1);
RED_ave = movmean(RED_raw, 5*FS, 1);

% hw_s = 10; % half-window width in second
% hw_sa = hw_s*FS; % half-window width in sample
% 
% % initialize the averaged arrays
% NIR_ave = zeros(length(NIR), 4);
% RED_ave = zeros(length(RED), 4);
% 
% for i =  - hw_sa:hw_sa
%     
%     NIR_ave = NIR_ave+circshift(NIR, i);
%     RED_ave = RED_ave+circshift(RED, i);
%     
% end
% 
% NIR_ave = NIR_ave./(2*hw_sa);
% RED_ave = RED_ave./(2*hw_sa);

% %% LOWPASS (BASELINE)
% FC = 0.1;
% 
% [b, a] = butter(4, FC/(FS/2), 'low'); % fc/(fs/2) is the cutoff freq in arc
% RED_f = filter(b, a, RED_raw_rev(:, PD_NUM));
% NIR_f = filter(b, a, NIR_raw_rev(:, PD_NUM));
% 
% % delete some head data point (since they yield wierd results)
% RED_f(1: START) = [];
% NIR_f(1: START) = [];
% TIME(length(TIME) - START + 1: end) = [];
% 
subplot(4, 2, 1);
plot(TIME, NIR_raw, 'g-',...
    TIME, NIR_ave, 'b--');
title('Filtered Data', 'fontsize', 30)
ylabel('Voltage(mV)', 'fontsize', 16)
set(gca,'FontSize', 14);
xlim([1, TIME(end)]);
% ylim([350, 420]);

subplot(4, 2, 3);
plot(TIME, RED_raw, 'm-',...
    TIME, RED_ave, 'r--');
xlabel('Time(s)', 'fontsize', 16)
ylabel('Voltage(mV)', 'fontsize', 16)
set(gca,'FontSize', 14);
xlim([1, TIME(end)]);
% ylim([380, 420]);

%% CALCULATE DIFFERENTIAL

% initial light intensity
RED_I0 = 390; % mV
NIR_I0 = 412; % mV

% absorption
RED_ab = log10(RED_I0./RED_ave);
NIR_ab = log10(NIR_I0./NIR_ave);

% differential
RED_ab_diff = circshift(RED_ab,-1)-RED_ab;
NIR_ab_diff = circshift(NIR_ab,-1)-NIR_ab;
%
subplot(4, 2, 2);
plot(TIME, NIR_ab, 'b',...
    TIME, 1000*NIR_ab_diff, 'g');
title('Absorption and differential', fontsize=30);
ylabel('D', fontsize=16);
set(gca, FontSize=14);
xlim([0, TIME(end)]);
ylim([-1e-3, 6e-3]);

subplot(4, 2, 4);
plot(TIME, RED_ab, 'r',...
    TIME, 1000*RED_ab_diff, 'm');
xlabel('Time(s)', fontsize=16);
ylabel('D', fontsize=16);
set(gca, FontSize=14);
xlim([0, TIME(end)]);
ylim([-1e-2, 4e-2]);

% subplot(4, 2, 2);
% plot(TIME, NIR_ab_diff, 'b');
% title('Absorption differential', fontsize=30);
% ylabel('D_ab (mm-1)', fontsize=16);
% set(gca, FontSize=14);
% xlim([0, TIME(end)]);
% % ylim([-5e-5, 5e-5]);
% 
% subplot(4, 2, 4);
% plot(TIME, RED_ab_diff, 'r');
% xlabel('Time(s)', fontsize=16);
% ylabel('D_ab (mm-1)', fontsize=16);
% set(gca, FontSize=14);
% xlim([0, TIME(end)]);
% % ylim([-5e-5, 5e-5]);

%% CALCULATE CONCENTRATION DIFFERENTIAL

extin_ox_RED = 0.011; % mm-1
extin_ox_NIR = 0.028; % mm-1
extin_deox_RED = 0.106; % mm-1
extin_deox_NIR = 0.018; % mm-1
dB = 2; % mm

conc_diff = (1/dB) .* mtimes(...
    inv([extin_deox_NIR, extin_ox_NIR;...
    extin_deox_RED, extin_ox_RED]),...
    [RED_ab_diff.'; NIR_ab_diff.']);
% shape of conc_diff: 2, length(TIME)
% conc_diff(1, :) stores deoxyhemoglobin conc change;
% conc_diff(2, :) stores oxyhemoglobin conc change

% subplot(4, 2, 2);
% plot(TIME, conc_diff(1, :), 'b')  % Raw Data Plot
% title('deox conc differential', fontsize=30)
% ylabel('mM', fontsize=16)
% set(gca, FontSize=14);
% xlim([0, TIME(end)]);
% ylim([-1e-5, 1e-5]);
% 
% subplot(4, 2, 4);
% plot(TIME, conc_diff(2, :), 'r')  % Raw Data Plot
% title('ox conc differential', fontsize=30)
% xlabel('Time(s)', fontsize=16)
% ylabel('mM', fontsize=16)
% set(gca, FontSize=14);
% xlim([0, TIME(end)]);
% ylim([-1e-5, 1e-5]);

%% CALCULATE OXYGEN SATURATION

% Initial concentration
conc_deox0 = 0.62; % mM
conc_ox0 = 1.86; % mM

conc_deox = conc_deox0;
conc_ox = conc_ox0;

tissue_SpO2 = zeros(length(TIME), 1);
conc = zeros(2, length(TIME));

for t = 1: length(TIME)
    tissue_SpO2(t) = 100*conc_ox/(conc_ox+conc_deox);
    conc(1, t) = conc_deox;
    conc(2, t) = conc_ox;
    
    conc_deox = conc_deox+conc_diff(1, t);
    conc_ox = conc_ox+conc_diff(2, t);
end

subplot(4, 2, 5);
plot(TIME, conc(1, :), 'b',...
    TIME, 1000*conc_diff(1, :), 'g');
title('deox conc and differential', fontsize=30);
ylabel('mM', fontsize=16);
set(gca, FontSize=14);
xlim([0, TIME(end)]);
ylim([-0.1, 0.3]);

subplot(4, 2, 7);
plot(TIME, conc(2, :), 'r',...
    TIME, 1000*conc_diff(2, :), 'm');
title('ox conc and differential', fontsize=30);
xlabel('Time(s)', fontsize=16);
ylabel('mM', fontsize=16);
set(gca, FontSize=14);
xlim([0, TIME(end)]);
ylim([-0.1, 1]);

subplot(2, 2, 4);
plot(TIME, tissue_SpO2, 'b')  % Raw Data Plot
xlabel('Time(s)', fontsize=16)
% ylabel('Voltage', 'fontsize', 16)
set(gca, FontSize=14);
xlim([0, TIME(end)]);
% ylim([85, 86]);

% subplot(4, 2, 8);
% plot(TIME, NIR_f, 'r')  % Raw Data Plot
% xlabel('Time(s)', 'fontsize', 16)
% % ylabel('Voltage', 'fontsize', 16)
% set(gca,'FontSize', 14);
% xlim([0, TIME(end)]);
% % ylim([80, 100]);

writematrix([TIME, tissue_SpO2], 'tissue_oximetry_res_hyp4_channel3.csv');
% save('tissue_oxygen.mat','tissue_SpO2');

