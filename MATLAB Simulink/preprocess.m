function [smooth_data, rect_data, data] = preprocess(data, fs, weight, window, F1, F2)
% This function is specific for LAB 3 
% With this function a signal can be preprocessed.
% The preprocess consists in bandpass filtering.
% Then, for the sake of visualization, there is a step of rectification and smoothing.

    % Bandpass filtering
    BP_filt = Cheb_38_filter(fs, F1, F2);
    data = filtfilthd(BP_filt, data);
    
   
    % FOR VISUALIZATION ONLY
    
    % TASK 2 - Rectify the signals
    % Where the signal has a negative value, modify to a correspondent positive value
    rect_data = abs(data);

    % Smooth the signals - Use a moving average filter with a fixed length window
    
    % Create the filter object
    movavgWindow = dsp.MovingAverage(window);
    % Apply the filter
    smooth_data = movavgWindow(rect_data);
    % Account for the delay introduced by FIR filtering
    smooth_data = cat(1,smooth_data(window/2+1:end),zeros(window/2, 1));
end

