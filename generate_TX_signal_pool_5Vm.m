% generate_TX_signals_pool.m
% Generates and organizes MOCZ signals for hardware transmission tests.
% We create 3 folders of packets ready for transmission - we will get 3 BER
% graphs as a result.
% -> packets with different K's
% -> packets with different L's
% -> packets with different SNR
clear all; close all; clc;

%% Setup Base Parameters and Directories
baseParams = makeParams(); % Load baseline parameters
baseParams.a = -5 ;
baseParams.b = 5 ;

num_signals_per_type = 5;  % Set how many signals to generate per variation

baseDir = fullfile(pwd, 'TX_signals/fiveVmax');
dirK   = fullfile(baseDir, '/MOCZ_packets/K_changes');
dirL   = fullfile(baseDir, '/MOCZ_packets/L_changes');
dirSNR = fullfile(baseDir, '/MOCZ_packets/SNR_changes');

% Create directories if they don't exist
if ~exist(baseDir, 'dir'), mkdir(baseDir); end
if ~exist(dirK, 'dir'),    mkdir(dirK);    end
if ~exist(dirL, 'dir'),    mkdir(dirL);    end
if ~exist(dirSNR, 'dir'),  mkdir(dirSNR);  end

fprintf('Starting TX Signal Generation Pool...\n');

%% Generate Signals for K Changes (5 to 20)
fprintf('\n--- Generating K variations ---\n');
K_vec = 5:20;
for k_val = K_vec
    for iter = 1:num_signals_per_type
        simParams = baseParams;
        simParams.K = k_val;
        simParams.power = 100;
        
        % Recalculate K-dependent variables!
        Kidxs = 0:k_val-1;
        simParams.R = sqrt(1 + 2 * simParams.lambda * sin(pi/k_val));
        simParams.theta_c = ((2*pi) * (Kidxs / k_val))';
        simParams.usedIdx = false(2^k_val, 1); % Reset tracking for each run
        
        % Generate
        [P, signal_pb_total, x_bb, group_delay, simParams] = generateFunc(simParams);
        
        % Save with iteration number
        fileName = sprintf('signal_Vm%d_K%d_L%d_Power100_run%d.mat', baseParams.b, k_val, baseParams.L, iter);
        save(fullfile(dirK, fileName), 'signal_pb_total', 'P', 'x_bb', 'simParams');
        fprintf('Saved: %s\n', fileName);
    end
end

rehash;   % Forces MATLAB to refresh its folder and path cache
pause(2); % Gives OneDrive 2 seconds to finish locking/syncing files

%% Generate Signals for L Changes (5ms to 50ms)
fprintf('\n--- Generating L variations ---\n');
L_vec = 1:1:50; % Guard time in symbols (which directly equals ms since Tsym = 1ms)
for l_val = L_vec
    for iter = 1:num_signals_per_type
        simParams = baseParams;
        simParams.L = l_val; 
        simParams.usedIdx = false(2^simParams.K, 1); % Reset tracking for each run
        simParams.power = 100;
        
        % Generate
        [P, signal_pb_total, x_bb, group_delay, simParams] = generateFunc(simParams);
        
        % Save with iteration number
        fileName = sprintf('signal_Vm%d_K%d_L%d_Power100_run%d.mat', baseParams.b, baseParams.K, l_val, iter);
        save(fullfile(dirL, fileName), 'signal_pb_total', 'P', 'x_bb', 'simParams');
        fprintf('Saved: %s (L = %d ms)\n', fileName, l_val);
    end
end

rehash;   % Forces MATLAB to refresh its folder and path cache
pause(2); % Gives OneDrive 2 seconds to finish locking/syncing files

%% Generate Signals for SNR Changes (Power Scaling)
fprintf('\n--- Generating SNR variations (Power Scaling) ---\n');
% Creates exactly 21 options between 0.5 and 1.0 (steps of 0.025)
power_scales = 0.5:0.025:1;
for p_scale = power_scales
    for iter = 1:num_signals_per_type
        simParams = baseParams;
        simParams.usedIdx = false(2^simParams.K, 1); % Reset tracking for each run
        
        % Generate standard signal
        [P, signal_pb_total, x_bb, group_delay, simParams] = generateFunc(simParams);
        
        % *** IMPORTANT MATH ***
        % To scale the POWER by p_scale, we must scale the AMPLITUDE by sqrt(p_scale)
        % Because power is proportional to the square of the amplitude
        amp_scale = sqrt(p_scale);
        x_bb_scaled = x_bb * amp_scale;
        signal_pb_total = signal_pb_total * amp_scale;
        
        % Multiply by 1000 for the filename to avoid decimal points in the string 
        % (e.g., 0.525 becomes 525)
        p_label = round(p_scale * 1000);
        simParams.power = p_label;
        
        % Save with iteration number
        fileName = sprintf('signal_Vm%d_K%d_L%d_Power%d_run%d.mat', baseParams.b, baseParams.K, baseParams.L, p_label, iter);
        save(fullfile(dirSNR, fileName), 'signal_pb_total', 'P', 'x_bb_scaled', 'simParams', 'p_scale', 'amp_scale');
        fprintf('Saved: %s (Amplitude Scale Factor: %.4f)\n', fileName, amp_scale);
    end
end

fprintf('\nDone! All signals generated and organized in TX_signals folder.\n');

%---------------------------------------------------------------------------------

%===========================================%
% Create signals with chirp at the begining
%===========================================%
%% Generate and Load the Chirp Signal
fprintf('Calling chirpFnc(1) to generate the chirp...\n');

%generate and save 'Chirp_Tx.mat'
chirpFnc(1); 

%ensure the file is written to the disk
pause(1); 

% Load the generated chirp file
if ~isfile('Chirp_Tx.mat')
    error('Chirp_Tx.mat was not found. Check if chirpFnc(1) executed correctly.');
end

chirpData = load('Chirp_Tx.mat');
tx_signal = chirpData.y_scaled; % chirpFnc(1) saves the signal as 'y_scaled'
Fs = chirpData.Fs;

% Optional: Add a tiny silence gap (e.g., 10ms) 
% between the chirp and MOCZ signal to prevent multipath overlap.
% gap_samples = zeros(0.01 * Fs, 1); 
% preamble = [tx_signal; gap_samples]; 
preamble = tx_signal;

%% Setup Directories
srcDir = fullfile(pwd, 'TX_signals');
destDir = fullfile(pwd, 'TX_signals_withChirp');

subfolders = {'K_changes', 'L_changes', 'SNR_changes'};

% Create the main destination folder
if ~exist(destDir, 'dir'), mkdir(destDir); end

%% Process Each Folder
for i = 1:length(subfolders)
    folderName = subfolders{i};
    srcSub = fullfile(srcDir, folderName);
    destSub = fullfile(destDir, folderName);
    
    % Create the subfolder in the new directory
    if ~exist(destSub, 'dir'), mkdir(destSub); end
    
    % Get all .mat files in the source subfolder
    files = dir(fullfile(srcSub, '*.mat'));
    
    fprintf('\n--- Processing folder: %s (%d files) ---\n', folderName, length(files));
    
    for j = 1:length(files)
        oldName = files(j).name;
        srcFile = fullfile(srcSub, oldName);
        
        % Load the original data into a struct
        data = load(srcFile);
        
        % Append the chirp preamble to the baseband signal.
        if isfield(data, 'x_bb_scaled')
            data.x_bb_scaled = [preamble; data.x_bb_scaled];
        elseif isfield(data, 'x_bb')
            data.x_bb = [preamble; data.x_bb];
        end
        
        % Also append to signal_pb_total to keep array lengths synchronized
        if isfield(data, 'signal_pb_total')
             data.signal_pb_total = [preamble; data.signal_pb_total];
        end
        
        % Create new filename with _Chirp appended
        [~, baseName, ext] = fileparts(oldName);
        newName = sprintf('%s_Chirp%s', baseName, ext);
        destFile = fullfile(destSub, newName);
        
        % Save the modified struct back to the new directory
        save(destFile, '-struct', 'data');
        fprintf('Saved: %s\n', newName);
    end
end

fprintf('\nDone! All chirped signals are ready in %s.\n', destDir);
