function [P_rec_D, P_rec_G, P_rec_M] = MOCZRealDecoding(simParams, SNR, x_pb, P)
% MOCZREALDECODING decodes a real recorded passband signal.
% It performs carrier down-conversion, matched filtering, exact
% Zadoff-Chu synchronization, downsampling, and DiZeT decoding.
% Note: 'P' is optional, used only for the plotting visualization.

    % -- Extract Parameters --
    K = simParams.K;
    B = simParams.B;
    L = simParams.L;
    Tsym = simParams.Tsym;
    Fs = simParams.Fs;
    fc = simParams.fc;
    beta = simParams.betha;
    sps = Tsym * Fs;
    
    P_rec_D = zeros(K, B);
    P_rec_G = zeros(K, B); % Initialize the second output for decoded messages
    P_rec_M = zeros(K, B);

    %% Adding channel
    %for real experiments remove this section. 
    h_channel = [0.8, -0.45, 0.2, -0.1 ,0.003]; % exmps: [0.8, -0.45, 0.2, -0.1 ,0.003], [0.8, -0.1, 0.003] 
    h_channel_upsmp = upsample(h_channel, sps);
    h_channel_upsmp = h_channel_upsmp / norm(h_channel_upsmp); %normalize

    %old:
    chanOutput = filter(h_channel_upsmp, 1, x_pb); %removes tail with last echos
    %new (chat suggested):
    % x_pb_padded = [x_pb; zeros(group_delay,1)];
    % chanOutput = filter(h_channel_upsmp, 1, x_pb_padded); %check it out!!
%     chanOutput = x_pb;
    
    x_pb_noisy = awgn(chanOutput, SNR ,'measured');

    x_pb = x_pb_noisy;

    %% Matched Filtering (Baseband Conversion at Fs)
    span = 6;
    h_coeff = rcosdesign(beta, span, sps);
    
    % Create time vector based on the actual length of the recording
    t = (0:length(x_pb)-1)' / Fs;
    
    % Downconvert to baseband
    mixed_x = x_pb .* exp(-1i*2*pi*fc*t);
    
    % Apply matched filter at High Rate (Fs)
    x_bb_filtered = filter(h_coeff, 1, mixed_x);

    %% Exact Zadoff-Chu Synchronization
    N = simParams.zadoffChuPair(2);
    u = simParams.zadoffChuPair(1);
    zc_ref = zadoffChuSeq(u, N);
    
    % Upsample reference to match Fs rate
    zc_upsamp = upsample(zc_ref, sps);

    % Cross-correlate filtered signal with the upsampled reference
    [xc, lags] = xcorr(x_bb_filtered, zc_upsamp);
    abs_xc = abs(xc);
    [maxVal, maxIdx] = max(abs_xc);
    peakLag = lags(maxIdx);
    
    % -- Peak Validation (PAPR check) --
    % If the peak isn't significantly higher than the average correlation,
    % you are likely just looking at noise.
    if (maxVal / mean(abs_xc)) < 15 
        warning('Low correlation peak detected. Packet may be missing or heavily corrupted.');
    end
    
    %%  Downsampling
    % The correlation peak implicitly accounts for the filter's group delay
    % peakLag tells us exactly where the first sample of the ZC sequence is.
    sync_start_sample_Fs = peakLag + 1;
    
    % To get to the first data symbol, we skip the ZC sequence (N) and the Guard (L)
    data_start_sample_Fs = sync_start_sample_Fs + ((N + L) * sps);
    
    % We need B messages, each of length (K + 1 + L)
    total_data_symbols = B * (K + 1 + L);
    end_sample_Fs = data_start_sample_Fs + (total_data_symbols - 1) * sps;
    
    % Safety check: Did the recording cut off before the packet ended?
    if end_sample_Fs > length(x_bb_filtered)
        error('Recorded signal is too short to contain the full B messages.');
    end
    
    % Downsample starting EXACTLY at the optimal phase
    sampled_data = x_bb_filtered(data_start_sample_Fs : sps : end_sample_Fs);
    x_decoded = sqrt(2) * sampled_data;
    
    %% 4. Decoding & Visualization
    x_decoded_matrix = reshape(x_decoded, K + 1 + L, B);
    
    R = simParams.R;
    theta_c = simParams.theta_c;

    for b = 1 : B
        x_vec_decoded = x_decoded_matrix(:, b);
        
        % Run decoder on the current message block
        message_rec = DiZeT(R, theta_c, x_vec_decoded, K);
        P_rec_D(:, b) = message_rec;
        message_rec = Greedy(x_vec_decoded, simParams);
        P_rec_G(:, b) = message_rec; % Store the greedy decoded message
        if K <= 8
            message_rec = ML(x_vec_decoded, simParams);
            P_rec_M(:, b) = message_rec;
        else
            fprintf('ML decoder was not activated because K > 8');
            P_rec_M(:, b) = zeros(K,1);
        end
        
        % Visualization (Only runs if figFlag is 1 AND 'P' was provided)
        if simParams.figFlag && nargin == 4
            decoded_alphas = roots(flip(x_vec_decoded.'));
            M = P(:, b);
            alphas = ((1-M)*(R^(-1)) + M * R) .* exp(1i*theta_c);

            figure;
            h1 = scatter(real(decoded_alphas), imag(decoded_alphas), 50, 'filled', 'b');
            grid on; axis equal; hold on;
            h2 = scatter(real(alphas), imag(alphas), 50, 'filled', 'g');
            
            thetas_for_plot = linspace(0, 2*pi, 300);
            plot(R * exp(1i*thetas_for_plot), 'r', 'LineWidth', 1);
            plot((1/R) * exp(1i*thetas_for_plot), 'r', 'LineWidth', 1);
            plot(exp(1i*thetas_for_plot), 'k--', 'LineWidth', 0.5);
            
            legend([h1, h2], {'Decoded Zeros', 'Original Zeros'});
            title(sprintf('Huffman BMOCZ Zeros - Message %d', b));
            xlabel('Re(z)'); ylabel('Im(z)');
            hold off;
        end
    end
end