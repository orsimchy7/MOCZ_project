
% analizing results of generated MOCZ signal - BER
clear all;
% close all;
%% load variables from mat file whose created by GenerateMOCZsignal
% load("L5B200K12N63u2.mat"); % 'signal_pb_total', 'P', 'simParams'
[simParams] = makeParams();

%% preparing
B = simParams.B;
K = simParams.K;
BW = 1 / simParams.Tsym; %assuming beta = 1
BWn = 3 * BW;
EbNo_dB = 0:18;
% SNRa = EbNo_dB + 10*log10(BW / BWn); %in BMOCZ, 1 symbol is 1 bit
% SNRa = -10 : 1 : 6;
SNRa = 10*log10((K/((K+1)*simParams.Tsym*BW*(1+simParams.betha))) * 10.^(EbNo_dB/10));
% SNRa = ones(19,1)*100;
% BER = zeros(length(SNRa), 1);
errorsNum_D = zeros(length(SNRa), 1);
% errorsNum_G = zeros(length(SNRa), 1);
errorsNum_M = zeros(length(SNRa), 1);
BER_D = zeros(length(SNRa), 1);
% BER_G = zeros(length(SNRa), 1);
BER_M = zeros(length(SNRa), 1);

attemptNum = 200 + 30 * (1:1:length(SNRa)); %200

%% running for different additive noise SNR
for i = 1 : length(SNRa)
    fprintf('i = %d \n', i);
    SNR = SNRa(i);
    for j = 1: attemptNum(i)
        [P, signal_pb_total, ~, ~, simParams] = generateFunc(simParams);
        simParams.P = P;
        %[P_rec] = MOCZsimChannelNdecoding(simParams, SNR, signal_pb_total, P, x_bb, group_delay);

        %New Decoder (includes zadoffchu crosscorelation)
        % [P_rec_D, P_rec_G, P_rec_M] = MOCZRealDecoding(simParams, SNR, signal_pb_total, P);
        % errorsNum_D(i) = errorsNum_D(i) + sum(abs(P - P_rec_D), 'all');
        % errorsNum_G(i) = errorsNum_G(i) + sum(abs(P - P_rec_G), 'all');
        % errorsNum_M(i) = errorsNum_M(i) + sum(abs(P - P_rec_M), 'all');
        % decode - old version only with dizet:
        % [P_rec_D, ~, ~] = decode_mocz(signal_pb_total, simParams, SNR);
        %visualize only last one:
        if j==attemptNum(i) && mod(i, 5) == 1
            simParams.figFlag = 1;
        end
        [P_rec_D, P_rec_ML, W_mat, ~, ~] = decode_mocz(signal_pb_total, simParams, SNR, i);
        errorsNum_D(i) = errorsNum_D(i) + sum(abs(P - P_rec_D), 'all');
        errorsNum_M(i) = errorsNum_M(i) + sum(abs(P - P_rec_ML), 'all');
        simParams.figFlag = 0;
    end
    BER_D(i) = errorsNum_D(i) / (B * K * attemptNum(i));
    % BER_G(i) = errorsNum_G(i) / (B * K * attemptNum(i));
    BER_M(i) = errorsNum_M(i) / (B * K * attemptNum(i));
    % fprintf('snr %d, BER_D is %d, BER_G is %d, BER_M is %d\n', SNR, BER_D(i), BER_G(i), BER_M(i));
    fprintf('snr %d, BER_D is %d, BER_M is %d\n', SNR, BER_D(i), BER_M(i));
    if mod(i, 5)== 1
        figure(2);
        j = 2 + (2/5)*(i - 1);
        subplot(4,2,j);
        imagesc(abs(W_mat)); %currently: always the last one
        colorbar
        title(sprintf('|(VV^H)^-^1| for uniform {0,1} zeros dist., R=%.2f, N/K=%.2f, EbN0=%.2f', simParams.R, simParams.N/K, EbNo_dB(i)))
        a=7;
    end
end

%% plotting BER results
% disp(errorsNum);
figure('Color', 'w'); % White background for reports
% Use semilogy for the logarithmic Y-axis
semilogy(EbNo_dB, BER_D, 'bo-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
hold on;
% semilogy(EbNo_dB(1:13), BER_G(1:13), 'ro-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
semilogy(EbNo_dB, BER_M, 'go-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'g');
grid on;
set(gca, 'YMinorGrid', 'on', 'XMinorGrid', 'off'); % Improves readability
xlabel('E_b/N_0 (dB)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Bit Error Rate (BER)', 'FontSize', 12, 'FontWeight', 'bold');
title('System BER Performance - normal channel and with zc', 'FontSize', 14);
% Optional: Set specific limits to make it look "standard"
% ylim([1e-6 1]); 
xlim([min(EbNo_dB) max(EbNo_dB)]);
% legend('decoder: DiZeT', 'decoder: Greedy', 'decoder: ML', 'Location', 'southwest');
legend('decoder: DiZeT', 'decoder: ML', 'Location', 'southwest');
hold off;
