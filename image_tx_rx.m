function [outputData, outParams] = image_tx_rx(mode, inputData, simParams)
    % IMAGE_TX_RX Handles MOCZ generation and decoding for images.
    %
    % TX MODE: outputData = signal to transmit | outParams = updated simParams
    % RX MODE: outputData = decoded 256x256 logical image matrix
    
    if strcmpi(mode, 'tx')
        % =========================================================
        % TRANSMIT MODE: Generate the signal
        % =========================================================
        % inputData is the 256x256 logical matrix from the GUI
        
        % Save matrix to a temp file so generateFunc can read it naturally
        temp_img_file = 'temp_gui_tx_image.png';
        imwrite(inputData, temp_img_file);
        
        if nargin < 3 || isempty(simParams)
            simParams = makeParams(); % Create default params if none provided
        end
        
        % Force image configuration parameters
        simParams.image_coding = 1;
        simParams.image_file = temp_img_file;
        simParams.K = 8; % K=8 is mandatory for 256x256
        simParams.a = 5; %top voltage amplitude [V]
        simParams.b = -5; %min voltage amplitude [V]
        
        % Dynamic Parameter Recalculation for K = 8
        K = simParams.K;
        L = simParams.L;
        N = K + L;
        simParams.N = N;
        
        Kidxs = 0:K-1;
        simParams.R = sqrt(1 + 2 * simParams.lambda * sin(pi/K));
        simParams.theta_c = ((2*pi) * (Kidxs / K))';
        simParams.usedIdx = false(2^K, 1);
        
        % Generate the actual transmission signal
        [P, signal_pb_total, x_bb, group_delay, simParams] = generateFunc(simParams);
        
        % Clean up the temporary file
        if isfile(temp_img_file)
            delete(temp_img_file);
        end
        
        outputData = signal_pb_total;
        outParams = simParams;
        
    elseif strcmpi(mode, 'rx')
        % =========================================================
        % RECEIVE MODE: Decode the hardware signal back into a photo
        % =========================================================
        x_pb = inputData(:); % Ensure input is a vertical column vector
        
        % -- Extract Parameters --
        K = simParams.K;
        B = simParams.B;
        L = simParams.L;
        Tsym = simParams.Tsym;
        Fs = simParams.Fs; 
        fc = simParams.fc;
        beta = simParams.betha;
        sps = Tsym * Fs;
        
        % ---------------------------------------------------------
        % Matched Filtering (Baseband Conversion)
        % ---------------------------------------------------------
        span = 6;
        h_coeff = rcosdesign(beta, span, sps);
        
        % Create time vector based on actual length of the recording
        t = (0:length(x_pb)-1)' / Fs;
        
        % Downconvert to baseband
        mixed_x = x_pb .* exp(-1i*2*pi*fc*t);
        
        % Apply matched filter at High Rate
        x_bb_filtered = filter(h_coeff, 1, mixed_x);
        
        % ---------------------------------------------------------
        % Exact Zadoff-Chu Synchronization
        % ---------------------------------------------------------
        N_zc = simParams.zadoffChuPair(2); % ZC Sequence Length
        u = simParams.zadoffChuPair(1);
        zc_ref = zadoffChuSeq(u, N_zc);
        
        % Upsample reference to match Fs rate and cross-correlate
        zc_upsamp = upsample(zc_ref, sps);
        [xc, lags] = xcorr(x_bb_filtered, zc_upsamp);
        abs_xc = abs(xc);
        [maxVal, maxIdx] = max(abs_xc);
        peakLag = lags(maxIdx);
        
        % Peak Validation (Warns you if the hardware didn't actually hear a signal!)
        if (maxVal / mean(abs_xc)) < 15 
            warning('Low correlation peak detected. Hardware recording may just be noise.');
        end
        
        % ---------------------------------------------------------
        % Downsampling
        % ---------------------------------------------------------
        sync_start_sample_Fs = peakLag + 1;
        data_start_sample_Fs = sync_start_sample_Fs + ((N_zc + L) * sps);
        
        % Calculate total length
        total_data_symbols = B * (K + 1 + L);
        end_sample_Fs = data_start_sample_Fs + (total_data_symbols - 1) * sps;
        
        if end_sample_Fs > length(x_bb_filtered)
            error('Hardware recording cut off early! Increase your listenTime padding in the GUI.');
        end
        
        % Downsample starting exactly at the optimal phase
        sampled_data = x_bb_filtered(data_start_sample_Fs : sps : end_sample_Fs);
        x_decoded = sqrt(2) * sampled_data;
        
        % Reshape into columns for the DiZeT loop
        x_decoded_matrix = reshape(x_decoded, K + 1 + L, B);
        
        % ---------------------------------------------------------
        % DiZeT Decoding Loop
        % ---------------------------------------------------------
        P_rec_D = zeros(K, B);
        R = simParams.R;
        theta_c = simParams.theta_c;
        
        for b = 1 : B
            x_vec_decoded = x_decoded_matrix(:, b);
            
            % Pass the chopped N-length block to your DiZeT function
            message_rec = DiZeT(R, theta_c, x_vec_decoded, K);
            P_rec_D(:, b) = message_rec;
        end
        
        % ---------------------------------------------------------
        % Reconstruct the Photo
        % ---------------------------------------------------------
        decoded_image = reshape(P_rec_D(:), [256, 256]);
        
        outputData = decoded_image;
        outParams = [];
        
    else
        error('Mode must be "tx" or "rx".');
    end
end