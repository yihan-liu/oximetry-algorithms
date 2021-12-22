% close all
clear;
clc;

cyc = 0;

%% READ IN
if cyc
    
    load('datasets/cycling_4.mat'); % load data
    
    TS = 1e-3; % interval for counting
    
    FS = 1/TS; % sampling frequency
    START = 5*FS;
    
    channel = 4;
    NIR_raw = NIR(:, channel);
    RED_raw = RED(:, channel);
    
    if length(NIR) > length(RED)
        
        NIR_raw(length(RED)+1:end) = [];
        TIME(length(RED)+1:end) = [];
        
    elseif length(NIR) < length(RED)
        
        RED_raw(length(NIR)+1:end) = [];
        TIME(length(NIR)+1:end) = [];
        
    end
    
    %% DENOISE (DC FILTER)
    
    hw_s = 20; % half-window width in second
    hw_sa = hw_s * FS; % half-window width in sample
    
    NIR_ave = movmean(NIR_raw, hw_sa);
    RED_ave = movmean(RED_raw, hw_sa);
    
    %% AC FILTER
    
    % bandpass filter
    heart_rate = [50 120];
    FC = heart_rate ./ 60;
    
    [b_ac, a_ac] = butter(2, FC/(FS/2), 'bandpass');
    RED_ac = filter(b_ac, a_ac, RED_raw);
    NIR_ac = filter(b_ac, a_ac, NIR_raw);
    
    % plot the dc data
    subplot(4, 2, 1);
    plot(TIME(ceil(START+1):end),...
        NIR_ave(ceil(START+1):end),...
        'b',...
        TIME(ceil(START+1):end),...
        NIR_ac(ceil(START+1):end)+93,...
        'c');
    title('Raw data & denoised data', 'fontsize', 30);
    ylabel('Voltage(mV)', 'fontsize', 16);
    set(gca,'FontSize', 14);
    
    subplot(4, 2, 3);
    plot(TIME(ceil(START+1):end),...
        RED_ave(ceil(START+1):end),...
        'r',...
        TIME(ceil(START+1):end),...
        RED_ac(ceil(START+1):end)+117,...
        'm');
    xlabel('Time(s)', 'fontsize', 16);
    ylabel('Voltage(mV)', 'fontsize', 16);
    set(gca,'FontSize', 14);
    
    %% RATIO & ENVELOPE
    
    RED_s = RED_ac ./ RED_ave;
    NIR_s = NIR_ac ./ NIR_ave;
    
    [NIR_upper, NIR_lower] = envelope(NIR_s, 0.8 * FS, 'peak');
    [RED_upper, RED_lower] = envelope(RED_s, 0.8 * FS, 'peak');
    
    NIR_amp = NIR_upper - NIR_lower;
    RED_amp = RED_upper - RED_lower;
    
    % plot the ratio
    subplot(4, 2, 2);
    plot(TIME(ceil(START+1):end),...
        NIR_amp(ceil(START+1):end),...
        'b',...
        TIME(ceil(START+1):end),...
        NIR_s(ceil(START+1):end),...
        'c',...
        TIME(ceil(START+1): end),...
        NIR_upper(ceil(START+1): end),...
        'g',...
        TIME(ceil(START+1): end),...
        NIR_lower(ceil(START+1): end),...
        'g');
    title('A_{NIR}(t)', 'fontsize', 30);
    set(gca,'FontSize', 14);
%     xlim([50 100])
    ylim([-5e-3, 5e-3]);
    
    subplot(4, 2, 4);
    plot(TIME(ceil(START+1): end),...
        RED_amp(ceil(START+1): end),...
        'r',...
        TIME(ceil(START+1): end),...
        RED_s(ceil(START+1): end),...
        'm',...
        TIME(ceil(START+1): end),...
        RED_upper(ceil(START+1): end),...
        'g',...
        TIME(ceil(START+1): end),...
        RED_lower(ceil(START+1): end),...
        'g');
    title('A_{RED}(t)', 'fontsize', 30);
    xlabel('Time(s)', 'fontsize', 16);
    set(gca,'FontSize', 14);
%     xlim([50 100])
    ylim([-5e-3, 5e-3]);
    
    %% CONVOLUTION
    
    log_ratio = log(NIR_amp + 1) ./ log(RED_amp + 1);
    
    width = 20*FS;
    integral_window = ones(width+1, 1);
    
    R = conv(log_ratio, integral_window, 'same') ./ width;
    
    subplot(2, 2, 3);
    plot(TIME(ceil(START+1):end),...
        R(ceil(START+1):end),...
        'b');
    title('R(t)', 'fontsize', 30);
    xlabel('Time(s)', 'fontsize', 16);
    set(gca,'FontSize', 14);
%     xlim([50 100]);
    
    %% BLOOD OXYGEN SATURATION
    
    extin_ox_RED = 0.011; % mm-1
    extin_ox_NIR = 0.028; % mm-1
    extin_deox_RED = 0.106; % mm-1
    extin_deox_NIR = 0.018; % mm-1
    
    % for cycling
    DPF_RED_NIR = 2.7; % mm
    
    % for hyp
    % DPF_RED_NIR = 1; % mm
    
    numerator = ...
        extin_deox_RED * DPF_RED_NIR * ones(length(TIME), 1) - ...
        extin_deox_NIR * R;
    denominator = ...
        (extin_deox_RED - extin_ox_RED) *...
        DPF_RED_NIR * ones(length(TIME), 1) + ...
        (extin_ox_NIR - extin_deox_NIR) * R;
    SpO2 = numerator ./ denominator;
    
    
    SpO2_ave = movmean(SpO2, 5*FS);
    
    subplot(2, 2, 4);
    plot(TIME(ceil(START+1):end),...
        SpO2_ave(ceil(START+1):end),...
        'b');
    xlabel('Time(s)', 'fontsize', 16)
    ylabel('SpO2(t)', 'fontsize', 16)
    title('Blood Oxygen Saturation', 'fontsize', 30)
    set(gca,'FontSize', 14);
    ylim([0.8, 1]);
%     xlim([50 100]);
    
else
    
    load('datasets/wrist_hypoxia_3.mat'); % load data
    
    TS = 1e-3; % interval for counting
    
    FS = 1/TS; % sampling frequency
    START = 15*FS;
    
    channel = 1;
    NIR_raw = NIR(:, channel);
    RED_raw = RED(:, channel);
    
    if length(NIR) > length(RED)
        
        NIR_raw(length(RED)+1:end) = [];
        TIME(length(RED)+1:end) = [];
        
    elseif length(NIR) < length(RED)
        
        RED_raw(length(NIR)+1:end) = [];
        TIME(length(NIR)+1:end) = [];
        
    end
    
    %% DENOISE
    
    hw_s = 0.1; % half-window width in second
    hw_sa = hw_s * FS; % half-window width in sample
    
    NIR_ave = movmean(NIR_raw, hw_sa);
    RED_ave = movmean(RED_raw, hw_sa);
    
    
    %% DENOISE (DC FILTER)
    
    hw_s = 1; % half-window width in second
    hw_sa = hw_s * FS; % half-window width in sample
    
    NIR_dc = movmean(NIR_ave, hw_sa);
    RED_dc = movmean(RED_ave, hw_sa);
    
    %% AC FILTER
    
    % bandpass filter
    heart_rate = [70 90];
    FC = heart_rate ./ 60;
    
    [b_ac, a_ac] = butter(2, FC/(FS/2), 'bandpass');
    RED_ac = filter(b_ac, a_ac, RED_raw);
    NIR_ac = filter(b_ac, a_ac, NIR_raw);
    
    % plot the dc data
    subplot(4, 2, 1);
    plot(TIME(ceil(START+1):end),...
        NIR_dc(ceil(START+1):end),...
        'b',...
        TIME(ceil(START+1):end),...
        NIR_ac(ceil(START+1):end)+NIR_dc(1)-NIR_ac(1),...
        'c',...
        TIME(ceil(START+1):end),...
        NIR_ave(ceil(START+1):end),...
        'g');
    title('Raw data & denoised data', 'fontsize', 30);
    ylabel('Voltage(mV)', 'fontsize', 16);
    set(gca,'FontSize', 14);
    
    subplot(4, 2, 3);
    plot(TIME(ceil(START+1):end),...
        RED_dc(ceil(START+1):end),...
        'r',...
        TIME(ceil(START+1):end),...
        RED_ac(ceil(START+1):end)+RED_dc(1)-RED_ac(1),...
        'm',...
        TIME(ceil(START+1):end),...
        RED_ave(ceil(START+1):end),...
        'g');
    xlabel('Time(s)', 'fontsize', 16);
    ylabel('Voltage(mV)', 'fontsize', 16);
    set(gca,'FontSize', 14);
    
    
    %% RATIO & ENVELOPE
    
    RED_s = RED_ac ./ RED_dc;
    NIR_s = NIR_ac ./ NIR_dc;
    
    [NIR_upper, NIR_lower] = envelope(NIR_s, 0.7 * FS, 'peak');
    [RED_upper, RED_lower] = envelope(RED_s, 0.7 * FS, 'peak');
    
    NIR_amp = NIR_upper - NIR_lower;
    RED_amp = RED_upper - RED_lower;
    
    % plot the ratio
    subplot(4, 2, 2);
    plot(TIME(ceil(START+1):end),...
        NIR_amp(ceil(START+1):end),...
        'b',...
        TIME(ceil(START+1):end),...
        NIR_s(ceil(START+1):end),...
        'c',...
        TIME(ceil(START+1): end),...
        NIR_upper(ceil(START+1): end),...
        'g',...
        TIME(ceil(START+1): end),...
        NIR_lower(ceil(START+1): end),...
        'g');
    title('A_{NIR}(t)', 'fontsize', 30);
    set(gca,'FontSize', 14);
    ylim([-1e-3, 1e-3]);
    
    subplot(4, 2, 4);
    plot(TIME(ceil(START+1): end),...
        RED_amp(ceil(START+1): end),...
        'r',...
        TIME(ceil(START+1): end),...
        RED_s(ceil(START+1): end),...
        'm',...
        TIME(ceil(START+1): end),...
        RED_upper(ceil(START+1): end),...
        'g',...
        TIME(ceil(START+1): end),...
        RED_lower(ceil(START+1): end),...
        'g');
    title('A_{RED}(t)', 'fontsize', 30);
    xlabel('Time(s)', 'fontsize', 16);
    set(gca,'FontSize', 14);
    ylim([-1e-3, 1e-3]);
    
    %% CONVOLUTION
    
    log_ratio = log(NIR_amp + 1) ./ log(RED_amp + 1);
    
    width = 15*FS;
    integral_window = ones(width+1, 1);
    
    R = conv(log_ratio, integral_window, 'same') ./ width;
    
    subplot(2, 2, 3);
    plot(TIME(ceil(START+1):end),...
        R(ceil(START+1):end),...
        'b');
    title('R(t)', 'fontsize', 30);
    xlabel('Time(s)', 'fontsize', 16);
    set(gca,'FontSize', 14);
    % xlim([100, 110]);
    
    %% BLOOD OXYGEN SATURATION
    
    extin_ox_RED = 0.011; % mm-1
    extin_ox_NIR = 0.028; % mm-1
    extin_deox_RED = 0.106; % mm-1
    extin_deox_NIR = 0.018; % mm-1
    
    DPF_RED_NIR = 1.4; % mm
    
    numerator = ...
        extin_deox_RED * DPF_RED_NIR * ones(length(TIME), 1) - ...
        extin_deox_NIR * R;
    denominator = ...
        (extin_deox_RED - extin_ox_RED) *...
        DPF_RED_NIR * ones(length(TIME), 1) + ...
        (extin_ox_NIR - extin_deox_NIR) * R;
    SpO2 = numerator ./ denominator;
    
    
    SpO2_ave = movmean(SpO2, 1*FS);
    
    subplot(2, 2, 4);
    plot(TIME(ceil(START+1):end),...
        SpO2_ave(ceil(START+1):end),...
        'b');
    xlabel('Time(s)', 'fontsize', 16)
    ylabel('SpO2(t)', 'fontsize', 16)
    title('Blood Oxygen Saturation', 'fontsize', 30)
    set(gca,'FontSize', 14);
%     ylim([0.8, 1.2]);
    
end
writematrix([TIME(ceil(START+1):end),...
    SpO2_ave(ceil(START+1):end)], 'SpO2_res_hyp3_ch2.csv');