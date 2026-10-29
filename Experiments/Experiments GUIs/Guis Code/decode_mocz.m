function [decoded_P, BER, x_decoded_matrix, W_mat] = decode_mocz(x_pb, simParams, dec_mode, n_trun)
% DECODE_HARDWARE_MOCZ Decodes real recorded MOCZ signals using only DiZeT.
% update: 1. dec_mode==1 is DiZet, dec_mode==2 is ML
% 2. added output W_mat [empty for dizet and not for ml].
        % ITS DIFFERENT FROM THE VERSION IN INTERNET COMPUTER!
        % since W_mat over here is a cell array in length B, not just the w_mat matrix
        % of last iteration. it contains all iterations here!
% 3. added n_trun that must satisfy: K+1 <= n_trun <= K+1+L TODOOOO
% Extracts the raw payload and calculates the exact Bit Error Rate.

    % -- Extract Parameters --
    K = simParams.K;
    B = simParams.B;
    L = simParams.L;
    Tsym = simParams.Tsym;
    Fs = simParams.Fs;
    fc = simParams.fc;
    beta = simParams.betha;
    sps = Tsym * Fs;
    
    x_pb = double(x_pb(:)); % Ensure it is a clean column vector
    
    % ---------------------------------------------------------
    % Matched Filtering (Baseband Conversion)
    % ---------------------------------------------------------
    span = 6;
    h_coeff = rcosdesign(beta, span, sps);
    t = (0:length(x_pb)-1)' / Fs;
    
    mixed_x = x_pb .* exp(-1i*2*pi*fc*t);
    x_bb_filtered = filter(h_coeff, 1, mixed_x);
    %no pulse shaping:
    % x_bb_filtered = mixed_x;

    % ---------------------------------------------------------
    % Exact Zadoff-Chu Synchronization
    % ---------------------------------------------------------
    N = simParams.zadoffChuPair(2);
    u = simParams.zadoffChuPair(1);
    zc_ref = zadoffChuSeq(u, N);
    
    zc_upsamp = upsample(zc_ref, sps);
    [xc, lags] = xcorr(x_bb_filtered, zc_upsamp);
    [maxVal, maxIdx] = max(abs(xc));
    peakLag = lags(maxIdx);
    
    if (maxVal / mean(abs(xc))) < 15 
        warning('Low correlation peak detected. Hardware recording may be noisy.');
    end
    
    % ---------------------------------------------------------
    % Downsampling
    % ---------------------------------------------------------
    sync_start_sample_Fs = peakLag + 1;
    data_start_sample_Fs = sync_start_sample_Fs + ((N + L) * sps);
    
    total_data_symbols = B * (K + 1 + L);
    end_sample_Fs = data_start_sample_Fs + (total_data_symbols - 1) * sps;
    
    if end_sample_Fs > length(x_bb_filtered)
        error('Recorded signal is too short to contain the full B messages.');
    end
    
    sampled_data = x_bb_filtered(data_start_sample_Fs : sps : end_sample_Fs);
    x_decoded = sqrt(2) * sampled_data;
    x_decoded_matrix = reshape(x_decoded, K + 1 + L, B);
    % update:
    if K+1 <= n_trun && n_trun <= K+1+L
        x_decoded_matrix = x_decoded_matrix(1: n_trun, :);
    end
    
    % ---------------------------------------------------------
    % DiZeT Decoding
    % ---------------------------------------------------------
    P_rec_D = zeros(K, B);
    R = simParams.R;
    theta_c = simParams.theta_c;
    % DO NOT DO THIS (BECUSE OF THE n_trun user update in GUI) - ASSUME IT APPLYED BEFORE!!!
    % [simParams] = runMLandGparams(simParams, K, L);
    W_mat = cell(B,1);
    for b = 1 : B
        x_vec_decoded = x_decoded_matrix(:, b);
        if dec_mode==1
            P_rec_D(:, b) = DiZeT(R, theta_c, x_vec_decoded, K);
            W_mat{b} = [];
        else % dec_mode==2
            [P_rec_D(:, b), W_mat{b}] = ML(x_vec_decoded, simParams);
        end
    end
    
    decoded_P = P_rec_D;
    
    % ---------------------------------------------------------
    % Bit Error Rate (BER) Calculation
    % ---------------------------------------------------------
    original_P = simParams.P; 
    
    % Compare the original matrix against the decoded matrix
    num_errors = sum(abs(original_P(:) - decoded_P(:)));
    total_bits = K * B;
    BER = num_errors / total_bits;
end