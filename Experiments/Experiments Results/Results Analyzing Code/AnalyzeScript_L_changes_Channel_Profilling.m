%% running L_changes
folderPath = '\\Recordings\fiveVmax\MOCZ_packets\L_changes';
files = dir(fullfile(folderPath, '*.mat'));

numRuns      = 5;
numL_values  = floor(length(files) / numRuns);
numDecoders  = 2; % 1 = DiZeT, 2 = ML
decoderNames = {'DiZeT', 'ML'};

L = nan(numL_values, 1);

% Running sums + per-L success counts. A skipped run is excluded from the
% average entirely (we divide by however many runs actually succeeded for
% that L), rather than silently averaging in an error as if it were data.
berL_sum      = zeros(numL_values, numDecoders);
berL_bit0_sum = zeros(numL_values, numDecoders);
berL_bit1_sum = zeros(numL_values, numDecoders);
runsUsed      = zeros(numL_values, 1);
skippedFiles  = {};

for k = 1 : numL_values
    for j = 1 : numRuns
        idx = numRuns * (k-1) + j;
        basefileName = files(idx).name;
        fullfileName = fullfile(folderPath, basefileName);
        fprintf('now: %s\n', basefileName);

        try
            load(fullfileName);

            if isnan(L(k))
                L(k) = currentParams.L;
            end

            DataMatrix  = [recordedData.Variables];
            signalData1 = DataMatrix(:, 1);
            signalData2 = DataMatrix(:, 2);

            % NOTE: decode_mocz.m now requires n_trun (K+1 <= n_trun <= K+1+L).
            % Using the full guard (K+1+L) preserves the old no-truncation
            % behavior AND is what makes Lt scale with L the way the bias
            % theory needs. Change this if you actually want to sweep n_trun
            % separately from L.
            n_trun = currentParams.K + 1 + currentParams.L;

            % Ground truth bits for this file (same message set for both
            % receivers). Assumes currentParams.P holds the transmitted bit
            % matrix, same convention as simParams.P elsewhere in the code --
            % double check this field exists in your saved recordings.
            trueBits = currentParams.P(:);
            mask0 = (trueBits == 0);
            mask1 = (trueBits == 1);

            for i = 1 : numDecoders
                [decoded_P1, error_rate1, ~] = decode_mocz(signalData1, currentParams, i, n_trun);
                [decoded_P2, error_rate2, ~] = decode_mocz(signalData2, currentParams, i, n_trun);

                errorR = (error_rate1 + error_rate2) * 0.5;
                berL_sum(k, i) = berL_sum(k, i) + errorR;

                % --- split by true bit value, pooling both receivers ---
                d1 = decoded_P1(:);
                d2 = decoded_P2(:);

                err0 = mean([abs(trueBits(mask0) - d1(mask0)); abs(trueBits(mask0) - d2(mask0))]);
                err1 = mean([abs(trueBits(mask1) - d1(mask1)); abs(trueBits(mask1) - d2(mask1))]);

                berL_bit0_sum(k, i) = berL_bit0_sum(k, i) + err0;
                berL_bit1_sum(k, i) = berL_bit1_sum(k, i) + err1;
            end

            runsUsed(k) = runsUsed(k) + 1;

        catch ME
            fprintf(2, '  -> SKIPPED %s  (%s)\n', basefileName, ME.message);
            skippedFiles{end+1} = basefileName; %#ok<SAGROW>
        end
    end
end

% Average using the actual number of successful runs per L (not a fixed
% numRuns). An L value with zero successful runs is left as NaN, which
% plot() renders as a gap in the line instead of a fake zero.
berL      = nan(numL_values, numDecoders);
berL_bit0 = nan(numL_values, numDecoders);
berL_bit1 = nan(numL_values, numDecoders);
for k = 1 : numL_values
    if runsUsed(k) > 0
        berL(k, :)      = berL_sum(k, :)      / runsUsed(k);
        berL_bit0(k, :) = berL_bit0_sum(k, :) / runsUsed(k);
        berL_bit1(k, :) = berL_bit1_sum(k, :) / runsUsed(k);
    else
        fprintf(2, 'WARNING: L = %d had ZERO successful runs out of %d -- plotted as a gap.\n', L(k), numRuns);
    end
end

% Files come back from dir() sorted alphabetically, so the L blocks were
% processed in an order like 1, 10, 11, ..., 19, 2, 20, ... instead of
% 1, 2, 3, ... Sort everything by the actual L value now so the graphs
% plot with a sensible, monotonic x-axis.
[L, sortIdx] = sort(L);
berL      = berL(sortIdx, :);
berL_bit0 = berL_bit0(sortIdx, :);
berL_bit1 = berL_bit1(sortIdx, :);
runsUsed  = runsUsed(sortIdx);

fprintf('\n=== Skip summary ===\n');
fprintf('%d of %d files skipped:\n', numel(skippedFiles), numel(files));
for s = 1 : numel(skippedFiles)
    fprintf('  - %s\n', skippedFiles{s});
end
fprintf('Runs used per L value:\n');
for k = 1 : numL_values
    fprintf('  L = %2d : %d / %d runs used\n', L(k), runsUsed(k), numRuns);
end
fprintf('\n');

%% Graph 1: BER split by true bit value, DiZeT vs ML side-by-side
% This is the primary evidence graph: watch for the bit=0 curve pulling
% away from bit=1 as L grows, and NOT closing back up.
colorZero = [0.85 0.20 0.20]; % red  -> bit = 0
colorOne  = [0.15 0.40 0.75]; % blue -> bit = 1
colorAgg  = [0.35 0.35 0.35]; % grey -> aggregate

figure('Name', 'BER split by true bit value', 'Color', 'w', 'Position', [100 100 1100 450]);
for i = 1 : numDecoders
    subplot(1, numDecoders, i);
    plot(L, berL_bit0(:, i)*100, '-o', 'LineWidth', 1.8, 'MarkerSize', 6, ...
        'Color', colorZero, 'MarkerFaceColor', colorZero);
    hold on;
    plot(L, berL_bit1(:, i)*100, '-s', 'LineWidth', 1.8, 'MarkerSize', 6, ...
        'Color', colorOne, 'MarkerFaceColor', colorOne);
    plot(L, berL(:, i)*100, '--^', 'LineWidth', 1.3, 'MarkerSize', 5, ...
        'Color', colorAgg, 'MarkerFaceColor', colorAgg);
    grid on; box on;
    xlabel('Number of guards L');
    ylabel('BER [%]');
    title(decoderNames{i});
    legend({'True bit = 0 (inside unit circle)', 'True bit = 1 (outside unit circle)', 'Aggregate'}, ...
        'Location', 'northwest');
    set(gca, 'FontSize', 11);
end
sgtitle('BER vs guard length L, split by transmitted bit value');
saveas(gcf, 'BER_vs_L_bitSplit.png');

%% Graph 2: the bias itself - gap between bit=0 and bit=1 BER, DiZeT vs ML
% Single clean number per L: if this line grows with L for DiZeT but
% stays flat (or grows much slower) for ML, that isolates the bias to
% DiZeT's decision rule rather than the zero-estimation/geometry itself.
gapDiZeT = (berL_bit0(:,1) - berL_bit1(:,1)) * 100;
gapML    = (berL_bit0(:,2) - berL_bit1(:,2)) * 100;

figure('Name', 'Decoder bias vs L', 'Color', 'w', 'Position', [100 100 700 450]);
plot(L, gapDiZeT, '-o', 'LineWidth', 2, 'MarkerSize', 6, 'Color', colorZero, 'MarkerFaceColor', colorZero);
hold on;
plot(L, gapML, '-s', 'LineWidth', 2, 'MarkerSize', 6, 'Color', colorOne, 'MarkerFaceColor', colorOne);
yline(0, '--k', 'LineWidth', 1);
grid on; box on;
xlabel('Number of guards L');
ylabel('BER(bit=0) - BER(bit=1)  [%]');
title('Zero-bit vs one-bit error asymmetry');
legend({'DiZeT', 'ML'}, 'Location', 'northwest');
set(gca, 'FontSize', 11);
saveas(gcf, 'BER_gap_vs_L.png');

%% Graph 3: aggregate BER, DiZeT vs ML
figure('Name', 'Aggregate BER vs L', 'Color', 'w', 'Position', [100 100 700 450]);
plot(L, berL(:,1)*100, '-o', 'LineWidth', 2, 'MarkerSize', 6, 'Color', colorZero, 'MarkerFaceColor', colorZero);
hold on;
plot(L, berL(:,2)*100, '-s', 'LineWidth', 2, 'MarkerSize', 6, 'Color', colorOne, 'MarkerFaceColor', colorOne);
grid on; box on;
xlabel('Number of guards L');
ylabel('BER [%]');
title('Aggregate BER vs guard length L');
legend(decoderNames, 'Location', 'northwest');
set(gca, 'FontSize', 11);
saveas(gcf, 'BER_aggregate_vs_L.png');

%% modeling the pool channel
load('\\wildwest\ecad\Matlab\projects\MOCZ_project_technion_orshoham\___all updated folders___\mocz_FilesFromPoolComputer\MOCZ_pool\NEW_Recordings_Pool_28_06\fiveVmax\Chirps\chirp_recorded_result.mat')
load('\\wildwest\ecad\Matlab\projects\MOCZ_project_technion_orshoham\___all updated folders___\mocz_FilesFromPoolComputer\MOCZ_pool\Chirp_Tx.mat');
refChirp = y_scaled;
ChirpsMatrix = [chirpRecordedData.Variables];
[corrChirp1, lags1] = xcorr(refChirp, ChirpsMatrix(:, 1));
[corrChirp2, lags2] = xcorr(refChirp, ChirpsMatrix(:, 2));

%% plotting chirps correlation
figure('Name', 'Chirp cross-correlation', 'Color', 'w', 'Position', [100 100 900 700]);
subplot(2,2,1);
plot(lags1, 10*log10(abs(corrChirp1)), 'LineWidth', 1.2, 'Color', colorZero);
grid on; box on; title('Rx1: |correlation| [dB]'); xlabel('lag'); ylabel('dB');

subplot(2,2,2);
plot(lags2, 10*log10(abs(corrChirp2)), 'LineWidth', 1.2, 'Color', colorOne);
grid on; box on; title('Rx2: |correlation| [dB]'); xlabel('lag'); ylabel('dB');

subplot(2,2,3);
plot(lags1, corrChirp1, 'LineWidth', 1.2, 'Color', colorZero);
grid on; box on; title('Rx1: correlation (linear)'); xlabel('lag');

subplot(2,2,4);
plot(lags2, corrChirp2, 'LineWidth', 1.2, 'Color', colorOne);
grid on; box on; title('Rx2: correlation (linear)'); xlabel('lag');

sgtitle('Chirp cross-correlation, both receivers');
saveas(gcf, 'Chirp_crossCorr.png');
