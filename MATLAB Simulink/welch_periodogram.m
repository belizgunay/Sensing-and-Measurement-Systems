function [Pall, freqs, MNF] = welch_periodogram(signal, fs, K, Nsegments, window, noverlap, NFFT)
% This function return the Welch periodogram

% INPUT ARGS
% signal: the signal to compute the periodogram
% fs: sampling frequency
% K: number of considered segments
% Nsegments: number of samples for each segment
% window: the window used in the welch approach
% noverlap: number of overlapping samples across consecutive segments
% NFFT: number of samples for the fft computation
% F1: inferior frequency limit
% F2: superior frequnecy limit

% OUTPUT ARGS
% Pall: matrix with PSDs of segments
% freqs: frequencies vector
% MNF: vector with MNF feature related to the observation instants

    for k = 1:K
        % Edge samples of the k-th segment
        K1 = (k-1)*Nsegments+1;
        K2 = k*Nsegments;

        % Define segment
        segment = signal(K1:K2);

        % Segment PSD
        [P, freqs] = pwelch(segment-mean(segment), window, noverlap, NFFT, fs);

        % TASK 5 - Extract the MNF feature (refer to slides)
        % Consider that freqs is the vector spanning from 0 to fs/2
        MNF(k) = sum(freqs.*P)/sum(P);

        % Save PSD of each segment
        Pall(k,:) = P;
    end
end

