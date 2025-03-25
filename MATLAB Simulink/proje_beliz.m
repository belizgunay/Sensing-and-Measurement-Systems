%% Preprocess the EMG Signal for Simulink Integration
clc;
clear;

% Sampling frequency
Fs = 2000; % [Hz]

%% Load Raw Data
load('C:\Users\beliz\Onedrive\Masaüstü\PROJECT_SAM\S1_DX_Inkjet_2Kg.mat'); % Ensure this file is uploaded
raw_signal = data; % Extract raw EMG signal

% Remove mean from the raw signal
raw_signal = raw_signal - mean(raw_signal);

% Time vector
t = (0:length(raw_signal)-1) / Fs; % Generate the time vector
t = t(:); % Ensure column vector
raw_signal = raw_signal(:); % Ensure column vector

%% Simulate SEN0240 Transfer Function
% Parameters of SEN0240
gain = 1000; % Gain of the sensor
offset_voltage = 1.5; % Offset voltage [V]
noise_std = 0.01; % Noise standard deviation (e.g., 10mV)

% Band-pass filter (15–250 Hz)
F1 = 15; % Lower cutoff frequency [Hz]
F2 = 250; % Upper cutoff frequency [Hz]
[b, a] = butter(4, [F1, F2] / (Fs / 2), 'bandpass'); % 4th-order Butterworth filter

% Apply the sensor transfer function
filtered_signal = filtfilt(b, a, raw_signal); % Filter the raw signal
sensor_output = gain * filtered_signal + offset_voltage; % Amplify and add offset
sensor_output = sensor_output + noise_std * randn(size(sensor_output)); % Add noise

% Save the processed signal for Simulink
processed_signal = sensor_output; % Match variable name expected by Simulink
save('processed_signal.mat', 'processed_signal', 't'); % Save to file
disp('Preprocessing complete. Processed signal saved as processed_signal.mat.');

%% Load Processed Signal for Simulink
load('processed_signal.mat'); % Ensure variable is in the workspace

%% Run the Simulink Model
simOut = sim('simulink_beliz.slx'); % Run your Simulink model and store the output in simOut

%% Continue Analysis...
% No changes here for further analysis steps

%% Extract Non-Buffered Output
non_buffered_data = simOut.emg_signal.Data; % Use 'Data' for timeseries structure
non_buffered_time = simOut.emg_signal.Time; % Use 'Time' for timeseries structure

%% Extract Buffered Output
buffered_data = simOut.simout.Data; % Use 'Data' for timeseries structure
buffered_time = simOut.simout.Time; % Use 'Time' for timeseries structure

% Ensure both outputs are column vectors
non_buffered_data = non_buffered_data(:);
buffered_data = buffered_data(:);

%% Analyze Non-Buffered Data
DUR = 1; % Segment duration [s]
Nsegments = DUR * Fs; % Samples per segment
K_non_buffered = floor(length(non_buffered_data) / Nsegments); % Number of segments (non-buffered)
welch_window = hamming(Nsegments);
noverlap = Fs / 3;
NFFT = 2^nextpow2(Nsegments);

%% Adjust Window for Welch's Method
segmentLength = Fs; % Use a segment length of 1 second
welch_window = hamming(segmentLength);
noverlap = segmentLength / 2; % 50% overlap
NFFT = 2^nextpow2(segmentLength);

%% Compute Welch's Periodogram for Non-Buffered Data
% Ensure enough data exists for at least one segment
if length(non_buffered_data) >= segmentLength
    [Pall_non, freqs] = pwelch(non_buffered_data, welch_window, noverlap, NFFT, Fs);
else
    error('Non-buffered data is too short for Welch analysis.');
end

%% Compute Welch's Periodogram for Buffered Data
% Reshape the buffered data into frames (if needed)
buffer_size = 64; % Match your buffer size in Simulink
num_frames = floor(length(buffered_data) / buffer_size); % Number of frames

if num_frames > 0
    buffered_frames = reshape(buffered_data(1:num_frames * buffer_size), buffer_size, num_frames);
    % Perform Welch's method on the concatenated buffered data
    [Pall_buf, ~] = pwelch(buffered_frames(:), welch_window, noverlap, NFFT, Fs);
else
    error('Buffered data is too short for Welch analysis.');
end

%% MNF and MDF Calculation for Non-Buffered Data
MNF_non = mean(freqs .* Pall_non, 'omitnan'); % Mean Frequency
cumsum_PSD_non = cumsum(Pall_non); % Cumulative PSD
MDF_non = freqs(find(cumsum_PSD_non >= sum(Pall_non) / 2, 1, 'first')); % Median Frequency

%% MNF and MDF Calculation for Buffered Data
MNF_buf = mean(freqs .* Pall_buf, 'omitnan'); % Mean Frequency
cumsum_PSD_buf = cumsum(Pall_buf); % Cumulative PSD
MDF_buf = freqs(find(cumsum_PSD_buf >= sum(Pall_buf) / 2, 1, 'first')); % Median Frequency

%% RMS Calculation
% Non-Buffered RMS
rms_non_buffered = rms(non_buffered_data);

% Buffered RMS (across frames)
rms_buffered = zeros(num_frames, 1);
for i = 1:num_frames
    frame = buffered_frames(:, i);
    rms_buffered(i) = rms(frame);
end

%% Display Results
fprintf('Feature Extraction Numerical Comparison:\n');
fprintf('----------------------------------------\n');
fprintf('RMS (Non-Buffered): %.4f\n', rms_non_buffered);
fprintf('RMS (Buffered): Mean = %.4f, Std Dev = %.4f\n', mean(rms_buffered), std(rms_buffered));
fprintf('----------------------------------------\n');
fprintf('MNF (Non-Buffered): %.4f\n', MNF_non);
fprintf('MNF (Buffered): %.4f\n', MNF_buf);
fprintf('----------------------------------------\n');
fprintf('MDF (Non-Buffered): %.4f\n', MDF_non);
fprintf('MDF (Buffered): %.4f\n', MDF_buf);

%% Validate Transfer Function: Plot Each Stage Separately

% Raw Signal
figure;
plot(t, raw_signal);
title('Raw EMG Signal (Input to SEN0240)');
ylabel('Amplitude (mV)');
xlabel('Time (s)');
grid on;

% Filtered Signal
figure;
plot(t, filtered_signal);
title('Filtered Signal (SEN0240: Band-Pass 15–250 Hz)');
ylabel('Amplitude (mV)');
xlabel('Time (s)');
grid on;

% Amplified Signal
amplified_signal = gain * filtered_signal; % Calculate amplified signal for visualization
figure;
plot(t, amplified_signal);
title('Amplified Signal (SEN0240 Transfer Function: Gain = 1000)');
ylabel('Amplitude (V)');
xlabel('Time (s)');
grid on;

% Offset Signal
offset_signal = amplified_signal + offset_voltage; % Add offset for visualization
figure;
plot(t, offset_signal);
title('Offset Signal (SEN0240 Transfer Function: Centered at 1.5 V)');
ylabel('Amplitude (V)');
xlabel('Time (s)');
grid on;

% Final Sensor Output (With Noise)
figure;
plot(t, sensor_output);
title('Final Sensor Output (SEN0240: With Noise)');
ylabel('Amplitude (V)');
xlabel('Time (s)');
grid on;

disp('Transfer function stages validated and plotted in separate figures.');

%% Filter the Non-Buffered and Buffered Signals (20–250 Hz)
% Apply the same band-pass filter to non-buffered and buffered signals
F1 = 20; % Lower cutoff frequency
F2 = 250; % Upper cutoff frequency
[b, a] = butter(4, [F1, F2] / (Fs / 2), 'bandpass'); % 4th-order Butterworth filter

% Filter the non-buffered signal
filtered_non_buffered = filtfilt(b, a, non_buffered_data);

% Filter the buffered signal
filtered_buffered = zeros(size(buffered_frames));
for i = 1:size(buffered_frames, 2)
    filtered_buffered(:, i) = filtfilt(b, a, buffered_frames(:, i));
end
filtered_buffered = filtered_buffered(:); % Flatten into a single column vector

%% Welch's Periodogram for Non-Buffered Signal
% Compute the PSD for the filtered non-buffered signal
[Pall_non, freqs] = pwelch(filtered_non_buffered, welch_window, noverlap, NFFT, Fs);

% Identify the valid range (non-zero PSD)
valid_idx_non = Pall_non > 0;
valid_PSD_non = Pall_non(valid_idx_non);
valid_freqs_non = freqs(valid_idx_non);

% Cumulative sum for valid PSD
cumsum_PSD_non = cumsum(valid_PSD_non);

% Compute MDF for non-buffered signal
MDF_non = valid_freqs_non(find(cumsum_PSD_non >= sum(valid_PSD_non) / 2, 1, 'first'));

% Compute MNF for non-buffered signal
MNF_non = mean(valid_freqs_non .* valid_PSD_non, 'omitnan');

%% Welch's Periodogram for Buffered Signal
% Compute the PSD for the filtered buffered signal
[Pall_buf, ~] = pwelch(filtered_buffered, welch_window, noverlap, NFFT, Fs);

% Identify the valid range (non-zero PSD)
valid_idx_buf = Pall_buf > 0;
valid_PSD_buf = Pall_buf(valid_idx_buf);
valid_freqs_buf = freqs(valid_idx_buf);

% Cumulative sum for valid PSD
cumsum_PSD_buf = cumsum(valid_PSD_buf);

% Compute MDF for buffered signal
MDF_buf = valid_freqs_buf(find(cumsum_PSD_buf >= sum(valid_PSD_buf) / 2, 1, 'first'));

% Compute MNF for buffered signal
MNF_buf = mean(valid_freqs_buf .* valid_PSD_buf, 'omitnan');

%% RMS Calculation
% Non-Buffered RMS
rms_non_buffered = rms(filtered_non_buffered);

% Buffered RMS (frame-by-frame)
rms_buffered = zeros(size(buffered_frames, 2), 1);
for i = 1:size(buffered_frames, 2)
    frame = filtered_buffered((i - 1) * size(buffered_frames, 1) + 1:i * size(buffered_frames, 1));
    rms_buffered(i) = rms(frame);
end

%% Display Results
fprintf('MDF and MNF Results:\n');
fprintf('----------------------------------------\n');
fprintf('MDF (Non-Buffered): %.4f Hz\n', MDF_non);
fprintf('MNF (Non-Buffered): %.4f Hz\n', MNF_non);
fprintf('MDF (Buffered): %.4f Hz\n', MDF_buf);
fprintf('MNF (Buffered): %.4f Hz\n', MNF_buf);

fprintf('----------------------------------------\n');
fprintf('RMS Results:\n');
fprintf('RMS (Non-Buffered): %.4f\n', rms_non_buffered);
fprintf('RMS (Buffered): Mean = %.4f, Std Dev = %.4f\n', mean(rms_buffered), std(rms_buffered));

%% Plot PSD for Comparison
% Non-Buffered PSD
figure;
plot(valid_freqs_non, valid_PSD_non, 'LineWidth', 1.5);
hold on;
yline(sum(valid_PSD_non) / 2, '--r', '50% Power', 'LineWidth', 1.2);
scatter(MDF_non, sum(valid_PSD_non) / 2, 80, 'r', 'filled');
title('Non-Buffered PSD (Filtered 20–250 Hz)');
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density');
legend('PSD', '50% Power', 'MDF');
xlim([0 100]); % Focus on 0–100 Hz
grid on;

% Buffered PSD
figure;
plot(valid_freqs_buf, valid_PSD_buf, 'LineWidth', 1.5);
hold on;
yline(sum(valid_PSD_buf) / 2, '--r', '50% Power', 'LineWidth', 1.2);
scatter(MDF_buf, sum(valid_PSD_buf) / 2, 80, 'r', 'filled');
title('Buffered PSD (Filtered 20–250 Hz)');
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density');
legend('PSD', '50% Power', 'MDF');
xlim([0 100]); % Focus on 0–100 Hz
grid on;

% Display the numerical MDF and MNF values
disp(['MDF (Non-Buffered Signal): ', num2str(MDF_non), ' Hz']);
disp(['MNF (Non-Buffered Signal): ', num2str(MNF_non), ' Hz']);
disp(['MDF (Buffered Signal): ', num2str(MDF_buf), ' Hz']);
disp(['MNF (Buffered Signal): ', num2str(MNF_buf), ' Hz']);
