function [decoded_P_D, decoded_P_ML, W_mat, BER, x_decoded_matrix] = decode_mocz(x_pb, simParams, snr, iter)
% DECODE_HARDWARE_MOCZ Decodes real recorded MOCZ signals using only DiZeT.
% Extracts the raw payload and calculates the exact Bit Error Rate for dizet.
% update: includes ML and returns random weight matrix (VV^H)^(-1).

    % -- Extract Parameters --
    K = simParams.K;
    B = simParams.B;
    L = simParams.L;
    Tsym = simParams.Tsym;
    Fs = simParams.Fs;
    fc = simParams.fc;
    beta = simParams.betha;
    sps = Tsym * Fs;
    
    % Adding channel
    %for real experiments remove this section. 

    % ----- constant channel------
    % h_channel = [0.75, -0.35, 0.1, -0.02, 0.003];
    % h_channel_upsmp = upsample(h_channel, sps);
    % h_channel_upsmp = h_channel_upsmp / norm(h_channel_upsmp); %normalize

    % ---------channel depended on L-------------
    % Create an exponentially decaying channel of length 2L:
    L_channel = 2 * simParams.L; 
    % Create an exponentially decaying channel of length L:
    L_channel = simParams.L;
    decay_factor = 0.5; % How fast the echoes fade - lower is faster decaying
    h_channel = decay_factor .^ (0:L_channel-1); 
    
    h_channel_upsmp = upsample(h_channel, sps);
    h_channel_upsmp = h_channel_upsmp / norm(h_channel_upsmp); % normalize

    chanOutput = filter(h_channel_upsmp, 1, x_pb); %removes tail with last echos
    %ber mode:
    SNR = snr;
    %test mode:
    % SNR = 20;%db
    x_pb_noisy = awgn(chanOutput, SNR ,'measured');

    x_pb = x_pb_noisy;
    
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
    data_start_sample_Fs = sync_start_sample_Fs + ((N + L-1) * sps);
    
    total_data_symbols = B * (K + L);
    end_sample_Fs = data_start_sample_Fs + (total_data_symbols - 1) * sps;
    
    if end_sample_Fs > length(x_bb_filtered)
        error('Recorded signal is too short to contain the full B messages.');
    end
    
    sampled_data = x_bb_filtered(data_start_sample_Fs : sps : end_sample_Fs);
    x_decoded = sqrt(2) * sampled_data;
    x_decoded_matrix = reshape(x_decoded, K + L, B);
    
    % ---------------------------------------------------------
    % DiZeT Decoding
    % ---------------------------------------------------------
    P_rec_D = zeros(K, B);
    R = simParams.R;
    theta_c = simParams.theta_c;

    for b = 1 : B
        x_vec_decoded = x_decoded_matrix(:, b);
        P_rec_D(:, b) = DiZeT(R, theta_c, x_vec_decoded, K);
        [P_rec_ML(:, b), G{b}] = ML(x_vec_decoded,simParams);
        vecs{b} = x_vec_decoded;
    end
    
    decoded_P_D = P_rec_D;
    decoded_P_ML = P_rec_ML;
    idx = randi(B);
    W_mat = G{idx};
    x_vec_decoded = vecs{idx}; % for visualize!

    
    % ---------------------------------------------------------
    % Bit Error Rate (BER) Calculation - only for DiZeT!
    % ---------------------------------------------------------
    original_P = simParams.P; 
    
    % Compare the original matrix against the decoded matrix
    num_errors = sum(abs(original_P(:) - decoded_P_D(:)));
    total_bits = K * B;
    BER = num_errors / total_bits;


    % Visualization (Only runs if figFlag is 1 AND 'P' was provided)
    figFlag = simParams.figFlag;
    if figFlag
        decoded_alphas = roots(flip(x_vec_decoded.'));
        P = simParams.P;
        M = P(:, b);
        alphas = ((1-M)*(R^(-1)) + M * R) .* exp(1i*theta_c);

        figure(2);
        j = 1 + (2/5)*(iter - 1);
        subplot(4,2,j);
        h1 = scatter(real(decoded_alphas), imag(decoded_alphas), 50, 'filled', 'b');
        grid on; axis equal; hold on;
        h2 = scatter(real(alphas), imag(alphas), 50, 'filled', 'g');
        
        thetas_for_plot = linspace(0, 2*pi, 300);
        plot(R * exp(1i*thetas_for_plot), 'r', 'LineWidth', 1);
        plot((1/R) * exp(1i*thetas_for_plot), 'r', 'LineWidth', 1);
        plot(exp(1i*thetas_for_plot), 'k--', 'LineWidth', 0.5);
        
        legend([h1, h2], {'Decoded Zeros', 'Original Zeros'});
        title(sprintf('Huffman BMOCZ Symbols. Message %d, snr=%.2f, K=%g, L=%g', idx, snr, K, L));
        xlabel('Re(z)'); ylabel('Im(z)');
        hold off;
    end
end
