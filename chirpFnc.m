function chirpFnc(RUN_MODE, recordedChirp_path)
% chirpFnc Runs the chirp generation, decoding, or simulation.
%
% Usage:
%   chirpFnc(1) - Generate Chirp (Save Tx .mat file)
%   chirpFnc(2) - Decode Chirp (Load recorded_chirp.mat & process)
%   chirpFnc(3) - Simulation (Pass through artificial channel & process)

    % Default to Simulation mode if no argument is provided
    if nargin < 1
        RUN_MODE = 3; 
    end

    % ==========================================
    % SYSTEM PARAMETERS
    % ==========================================
    T = 50e-3;          % duration of chirp [s]
    Fs = 200e3;         % sample frequency [Hz]
    t = (1/Fs:1/Fs: T)';
    f0 = 20000; %0.05*Fs;       % start frequency (10 kHz)
    f1 = 45000; %0.4*Fs;        % end frequency (80 kHz)
    V_peak = 10;        % Stretch for transmitter range [20Vpp]

    % Generate the signal
    y_ref = chirp(t, f0, T, f1, 'logarithmic');
    tx_signal = y_ref * V_peak;


    % ==========================================
    % MODE 1: GENERATE & SAVE ONLY
    % ==========================================
    if RUN_MODE == 1
        disp('--- MODE 1: GENERATING CHIRP ---');
        filename = 'Chirp_Tx.mat';
        y_scaled = tx_signal; 
        save(filename, 'y_scaled', 'Fs');
        disp(['Success: Chirp signal saved as ', filename]);
        
        % Verify visually
        figure;
        pspectrum(y_scaled, Fs, 'spectrogram', 'TimeResolution', 1e-3, ...
            'OverlapPercent', 99, 'Leakage', 0.85);
        title('Transmitted Signal Spectrogram');
        
        return; % Exits the function early
    end


    % ==========================================
    % MODE 2 & 3: SIGNAL SETUP
    % ==========================================
    has_true_taps = false; % Flag for plotting

    if RUN_MODE == 2
        disp('--- MODE 2: DECODING REAL CHANNEL ---');

        if nargin<2 || isempty(recordedChirp_path)
            disp('recorded chirp file not found. searching for default file');
            recordedChirp_path = 'chirp_recorded_result.mat';
        end
        
        % Load the recorded received signal
        if ~isfile(recordedChirp_path)
            error('%s not found in current directory.', recordedChirp_path);
        end
        rx_data = load(recordedChirp_path);
        
        % Dynamically grab the first variable inside recorded_chirp.mat 
        fields = fieldnames(rx_data);
        recorded_table = rx_data.(fields{1}); 

        if istable(recorded_table) || istimetable(recorded_table)
            if ismember('Dev5_ai0', recorded_table.Properties.VariableNames)
                rx_signal = recorded_table.Dev5_ai0;
                disp('Extracted data from channel Dev5_ai0');
            else
                error('Column: Dev5_ai0 was not found in the data table');
            end
        else
            rx_signal = recorded_table;
        end

        rx_signal = double(rx_signal(:));

    elseif RUN_MODE == 3
        disp('--- MODE 3: SIMULATION ---');
        
        %Create the 5-Tap Decaying Channel
        tap_delays = [0, 50, 120, 200, 300]; 
        tap_delays_samples = round(tap_delays * 1e-6 * Fs);
        tap_amps = [1.0, 0.6, 0.3, 0.15, 0.05];
        
        channel_length = max(tap_delays_samples) + 1;
        h_true = zeros(channel_length, 1);
        h_true(tap_delays_samples + 1) = tap_amps; 
        
        % 2. Pass signal through channel and add noise
        rx_signal = filter(h_true, 1, tx_signal);
        rx_signal = awgn(rx_signal, 20, 'measured'); 
        
        has_true_taps = true; % Tell the plotting section to draw the red stems
    end


    % ==========================================
    % UNIFIED DECODE & PROCESSING (Runs for Mode 2 & 3)
    % ==========================================

    %Matched Filter (Cross-Correlation)
    [corr_out, lags] = xcorr(rx_signal, tx_signal); %for decode mode tx_signal = the original clean chirp

    % Convert lags to microseconds and keep only positive delays
    time_lags_us = (lags / Fs) * 1e6;
    positive_idx = find(lags >= 0);
    corr_positive = corr_out(positive_idx);
    time_positive_us = time_lags_us(positive_idx);

    % Visualization: Impulse Response
    figure;
    plot(time_positive_us, abs(corr_positive), 'LineWidth', 1.5);
    hold on;

    if has_true_taps
        % Overlay true taps for simulation mode
        stem(tap_delays, tap_amps * max(abs(corr_positive)), 'r', 'LineWidth', 1.5);
        legend('Matched Filter Output', 'True Channel Taps');
    else
        legend('Matched Filter Output');
    end

    xlim([0 400]); 
    title('Estimated Channel Impulse Response (Matched Filter)');
    xlabel('Delay (\mu s)');
    ylabel('Correlation Magnitude');
    grid on;

    %Extract Peaks (Taps)
    threshold = 0.1 * max(abs(corr_positive)); 
    min_dist_samples = round(14e-6 * Fs); % Based on 70kHz BW resolution

    [peak_amps, peak_locs] = findpeaks(abs(corr_positive), ...
        'MinPeakHeight', threshold, 'MinPeakDistance', min_dist_samples);

    % Normalize delays (n=0 for first arrival)
    delay_samples = peak_locs - peak_locs(1);

    % Build discrete polynomial
    max_delay = max(delay_samples);
    h_discrete = zeros(max_delay + 1, 1);
    h_discrete(delay_samples + 1) = peak_amps; 

    %Visualization: Z-Plane
    figure;
    % Transpose h_discrete to a row vector to find roots correctly
    zplane(h_discrete', 1); 
    title('Z-Plane: Zeros of the Estimated Channel');

end