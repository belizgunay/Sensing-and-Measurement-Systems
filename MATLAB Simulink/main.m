%% LAB 2 
clc
clear all
close all

%% Load dataset

% Sampling frequency
fs = 2000; % [Hz]

% Select load
weight = 8;

path = '../Data/S1_inkjet/S1_DX_Inkjet_' + string(weight) + 'Kg.mat';
load_struct = load(path);

raw_data = load_struct.data;

clear load_struct, clear path

% Remove mean from data
raw_data = raw_data - mean(raw_data);

% Time vector
t = (1/fs:1/fs:length(raw_data)/fs)';

figure
subplot(1,2,1); plot(t, raw_data)
xlabel('Time [s]'), ylabel('Amplitude [mV]')
title("Time Domain")

fft_data = 1/fs*fft(raw_data);
NFFT = length(fft_data);
F = fs/NFFT;

% Frequencies vector
f = F*(0:NFFT-1);

subplot(1,2,2); plot(f, abs(fft_data).^2)
xlabel('Frequency [Hz]'), ylabel('Power [mV^2/Hz]')
xlim([0 fs/2])
title("Frequency Domain")

%suptitle('sEMG (Zero mean)')
sgtitle('sEMG (Zero mean)')

%% TASK 1 - Set the frequency limits [Hz]
F1 = 15;
F2 = 250;

%% TASK 2 - Inside function "preprocess"
% You can leave this space empty. You have to complete the function preprocess.m that is called in the next task

%% TASK 3 - Define the window for the moving average filter
len_window = fs/2; % samples

[smooth_data, rect_data, data] = preprocess(raw_data, fs, weight, len_window, F1, F2);

% Time vector
%t = (1/fs:1/fs:length(smooth_data)/fs)';

figure
plot(t, rect_data)
grid, hold on
plot(t, smooth_data)
xlabel('Time [s]'), ylabel('Amplitude [mV]')
title("sEMG - Time Domain")
legend('Rectified', 'Rectified and smoothed')

%% Limit the analysis to the interval with the contraction (30s by protocol design)

T = 30; % [s]

% TASK 4 - From the previous plot identify the instants of contraction start/stop 
T1 = 3;  % [s]

T2 = T1 + T; % [s]

t_contraction = t(T1*fs+1:T2*fs); 
signal = data(T1*fs+1:T2*fs);

% Segmentation in 1 second segments
DUR = 1; % [s]
Nsegments = DUR*fs; % samples
K = T*fs/Nsegments; % trials in the entire signal

%% Welch Periodogram 

% Set Welch's periodogram parameters
welch_window = hamming(Nsegments/DUR);
noverlap = fs/3;
NFFT     = 2.^nextpow2(Nsegments/DUR);

% Compute the periodogram
% TASK 5 - Inside function "welch_periodogram"
[Pall,freqs, MNF] = welch_periodogram(signal, fs, K, Nsegments, welch_window, noverlap, NFFT);

% Plot avg PSD
figure
area(freqs, mean(Pall), 'linewidth', 1)
grid
xlim([F1 F2])
xlabel('Frequency [Hz]'), ylabel('Power [mV^2/Hz]')
title('sEMG - Frequency Domain')

%% Assess fatigue
c     = polyfit(1:K, MNF, 1);
MNF_est = polyval(c,1:K);
clear c

figure
plot((1:K)*2, MNF, 'ks-')
grid, hold on
plot((1:K)*2, MNF_est, 'r--')
xlabel('Observation times [s]'), ylabel('MNF [Hz]')
title('Fatigue assessment - ' + string(weight) + ' Kg')

% TASK 6 - Describe the previous plot in your report
% Which is the behaviour of the MNF feature? Did you expect it?

%% COMPARE WITH ANOTHER ACQUISITION SESSION
%  The experimental protocol is assumed to be the same
new_weight = 2;

path = '../Data/S1_inkjet/S1_DX_Inkjet_' + string(new_weight) + 'Kg.mat';
load_struct = load(path);

new_raw_data = load_struct.data;

clear load_struct, clear path

% Remove mean from data
new_raw_data = new_raw_data - mean(new_raw_data);

[new_smooth_data, new_rect_data, new_data] = preprocess(new_raw_data, fs, new_weight, len_window, F1, F2);

% Time vector
new_t = 1/fs:1/fs:length(new_smooth_data)/fs;

figure
plot(new_t, new_rect_data)
grid, hold on
plot(new_t, new_smooth_data)
xlabel('Time'), ylabel('Amplitude [mV]')
title('sEMG ' + string(new_weight) + ' Kg - Time Domain')
legend('Rectified', 'Rectified and smoothed')


%% Contraction level comparison

figure
plot(t_contraction, smooth_data(T1*fs+1:T2*fs), 'linewidth', 2)
grid, hold on
plot(t_contraction, new_smooth_data(T1*fs+1:T2*fs), 'linewidth', 2)
xlabel('Time'), ylabel('Amplitude [mV]')
xlim([T1 T2])
title('Comparison sEMGs - Time Domain')
legend(string(weight) + ' Kg', string(new_weight)+ ' Kg')

% TASK 7 - Describe the previous plot in your report
% Which are your conclusions from the contraction level comparison? What did you expect?

%% Avg PSD comparison

new_signal = new_data(T1*fs+1:T2*fs);

% Compute the periodogram
[new_Pall,new_freqs, new_MNF] = welch_periodogram(new_signal, fs, K, Nsegments, welch_window, noverlap, NFFT);

figure
area(freqs, mean(Pall), 'linewidth', 1)
grid, hold on
area(new_freqs, mean(new_Pall), 'linewidth', 1)
xlim([F1 F2])
xlabel('Frequency [Hz]'), ylabel('Power [mV^2/Hz]')
title('Comparison sEMGs - Frequency Domain')
legend(string(weight) + ' Kg', string(new_weight)+ ' Kg')

% TASK 8 - Describe the previous plot in your report
% Which are your conclusions from the average PSD comparison? What did you expect?

%% Fatigue comparison
new_c     = polyfit(1:K, new_MNF, 1);
new_MNF_est = polyval(new_c,1:K);
clear new_c

figure
plot((1:K)*2, MNF, 'bs-', 'HandleVisibility','off')
grid, hold on
plot((1:K)*2, new_MNF, 'rs-', 'HandleVisibility','off')
plot((1:K)*2, MNF_est, 'b--', 'DisplayName', string(weight) + ' Kg')
plot((1:K)*2, new_MNF_est, 'r--', 'DisplayName', string(new_weight) + ' Kg')
legend show
xlabel('Observation times [s]'), ylabel('MNF [Hz]')
title('Fatigue comparison')

% TASK 9 - Describe the previous plot in your report
% Which are your conclusions from the fatigue comparison? What did you expect?
