%% SNR_changes analysis: BER vs SNR, at several EFFECTIVE guard lengths
%
% All recordings in this folder were made with the SAME physical guard
% length (currentParams.L -- currently 25 for all of them). Instead of
% needing separate recordings per L, we reuse each recording and re-decode
% it multiple times with different n_trun values, truncating how many
% guard samples the decoder actually gets to use:
%     n_trun = K + 1 + targetL
% "effective L" (targetL) must be <= the physically recorded L -- you
% can't truncate into existence guard samples that were never recorded.
%
% Same defensive logic as before: every file's processing is wrapped in
% try/catch, and it's committed to the running totals only if the WHOLE
% file (every target L, every decoder, both receivers) succeeds -- so a
% mid-file failure can't leave some target-L bins with one more sample
% than others.
%
% TODO: point this at wherever your SNR_changes recordings actually live.
folderPath = '\\wildwest\ecad\Matlab\projects\MOCZ_project_technion_orshoham\___all updated folders___\mocz_FilesFromPoolComputer\MOCZ_pool\NEW_Recordings_Pool_28_06\fiveVmax\Difficult_channel\MOCZ_packets\SNR_changes';
files = dir(fullfile(folderPath, '*.mat'));

numDecoders  = 2; % 1 = DiZeT, 2 = ML
decoderNames = {'DiZeT', 'ML'};

% Effective (truncated) guard lengths to test -- NOT the physical L the
% packets were recorded with, just how many guard samples n_trun asks the
% decoder to actually consume. All plotted together on one axes now.
Lvalues_toPlot = [1 4 15];
allLvalues     = sort(unique(Lvalues_toPlot));

%% Pass 1: load once per file, decode once per (targetL, decoder), accumulate
% key = 'L<target>_P<powerLabel>' -> struct('L',.., 'powerLabel',..,
%                                            'sumBER', 1xnumDecoders, 'count', scalar)
resultsMap     = containers.Map();
skippedFiles   = {};
warnedInvalidL = [];

% ML's decoder matrices (V_all/Cinv_all/combinations) are sized for a
% fixed N = K + L and depend only on K, R, theta_c -- not on any specific
% recording. Precompute them ONCE per requested effective L (not once per
% file) and reuse. DiZeT needs none of this; it's built fresh from the
% truncated vector every call regardless of L.
mlFieldsPerL = containers.Map('KeyType', 'double', 'ValueType', 'any');

for idx = 1 : numel(files)
    basefileName = files(idx).name;
    fullfileName = fullfile(folderPath, basefileName);
    fprintf('now: %s\n', basefileName);

    try
        load(fullfileName);

        actualL    = currentParams.L;   % physically recorded guard length
        powerLabel = currentParams.power;

        if mlFieldsPerL.Count == 0
            % One-time precompute, using this (or any) file's K/R/theta_c.
            for tL = allLvalues
                mlParams = runMLandGparams(currentParams, currentParams.K, tL);
                mlFieldsPerL(tL) = struct('N', mlParams.N, 'V_bank', mlParams.V_bank, ...
                    'combinations', mlParams.combinations, 'Cinv_all', mlParams.Cinv_all, ...
                    'V_all', mlParams.V_all);
            end
        end

        DataMatrix  = [recordedData.Variables];
        signalData1 = DataMatrix(:, 1);
        signalData2 = DataMatrix(:, 2);

        validTargets   = allLvalues(allLvalues <= actualL);
        invalidTargets = setdiff(allLvalues, validTargets);
        newlyInvalid   = setdiff(invalidTargets, warnedInvalidL);
        if ~isempty(newlyInvalid)
            fprintf(2, 'NOTE: requested L = %s exceed the physically recorded L = %d -- skipping those.\n', ...
                mat2str(newlyInvalid), actualL);
            warnedInvalidL = union(warnedInvalidL, newlyInvalid);
        end

        fileUpdates = containers.Map(); % this file's contribution -- merged only if the whole file succeeds

        for targetL = validTargets
            n_trun = currentParams.K + 1 + targetL;
            berThisTargetL = zeros(1, numDecoders);

            % Same physical L (needed for decode_mocz's sync/reshape), but
            % with ML's matrices swapped for ones sized to this targetL.
            mlFields = mlFieldsPerL(targetL);
            paramsForThisL = currentParams;
            paramsForThisL.N            = mlFields.N;
            paramsForThisL.V_bank       = mlFields.V_bank;
            paramsForThisL.combinations = mlFields.combinations;
            paramsForThisL.Cinv_all     = mlFields.Cinv_all;
            paramsForThisL.V_all        = mlFields.V_all;

            for i = 1 : numDecoders
                [~, error_rate1, ~] = decode_mocz(signalData1, paramsForThisL, i, n_trun);
                [~, error_rate2, ~] = decode_mocz(signalData2, paramsForThisL, i, n_trun);
                berThisTargetL(i) = (error_rate1 + error_rate2) * 0.5;
            end

            key = sprintf('L%d_P%d', targetL, powerLabel);
            fileUpdates(key) = struct('L', targetL, 'powerLabel', powerLabel, 'ber', berThisTargetL);
        end

        % Whole file succeeded for every requested target L -- now merge
        % its contribution into the running totals.
        updKeys = keys(fileUpdates);
        for kk = 1 : numel(updKeys)
            upd = fileUpdates(updKeys{kk});
            if isKey(resultsMap, updKeys{kk})
                entry = resultsMap(updKeys{kk});
            else
                entry = struct('L', upd.L, 'powerLabel', upd.powerLabel, ...
                                'sumBER', zeros(1, numDecoders), 'count', 0);
            end
            entry.sumBER = entry.sumBER + upd.ber;
            entry.count  = entry.count + 1;
            resultsMap(updKeys{kk}) = entry;
        end

    catch ME
        fprintf(2, '  -> SKIPPED %s  (%s)\n', basefileName, ME.message);
        skippedFiles{end+1} = basefileName; %#ok<SAGROW>
    end
end

fprintf('\n=== Skip summary ===\n');
fprintf('%d of %d files skipped:\n', numel(skippedFiles), numel(files));
for s = 1 : numel(skippedFiles)
    fprintf('  - %s\n', skippedFiles{s});
end
fprintf('\n');

%% Unpack the map into flat arrays and average
allValues = values(resultsMap);

Lall     = cellfun(@(v) v.L, allValues);           % effective (truncated) L
Pall     = cellfun(@(v) v.powerLabel, allValues);
countAll = cellfun(@(v) v.count, allValues);
berAll   = nan(numel(allValues), numDecoders);
for k = 1 : numel(allValues)
    if allValues{k}.count > 0
        berAll(k, :) = allValues{k}.sumBER / allValues{k}.count;
    end
end

% Relative SNR in dB -- see note in the previous version of this script:
% it's relative to whatever the ambient noise floor was during recording,
% not an absolute SNR.
snrDbAll = 10 * log10(Pall / 1000);

fprintf('Runs per (effective L, SNR) point:\n');
uniqCombos = sortrows(unique([Lall(:), Pall(:)], 'rows'), [1 2]);
for k = 1 : size(uniqCombos, 1)
    thisL = uniqCombos(k,1); thisP = uniqCombos(k,2);
    thisCount = countAll(Lall == thisL & Pall == thisP);
    fprintf('  L = %2d, Power = %4d : %d runs used\n', thisL, thisP, thisCount);
end
fprintf('\n');

%% Debug: print the actual BER numbers, per L, per SNR, per decoder.
% Use this to confirm whether a "missing" line in a plot is really
% missing data, vs. just plotted-at-zero-and-covered-by-a-later-line.
fprintf('=== BER by effective L, sorted by SNR ===\n');
for Lval = allLvalues
    rows = (Lall == Lval);
    if ~any(rows)
        fprintf('  L = %2d : NO DATA AT ALL (never made it into resultsMap)\n', Lval);
        continue;
    end
    [snrSorted, order] = sort(snrDbAll(rows));
    berRows = berAll(rows, :);
    berRows = berRows(order, :);
    fprintf('  L = %2d :\n', Lval);
    for kk = 1 : numel(snrSorted)
        fprintf('      SNR = %6.2f dB   DiZeT = %6.3f%%   ML = %6.3f%%\n', ...
            snrSorted(kk), berRows(kk,1)*100, berRows(kk,2)*100);
    end
end
fprintf('\n');

%% Color scale shared across BOTH panels and BOTH figures: bigger L = further along the map
cmap = turbo(numel(allLvalues)); % swap for parula(...) if turbo isn't available in your MATLAB version
colorForL = containers.Map('KeyType', 'double', 'ValueType', 'any');
for idx = 1 : numel(allLvalues)
    colorForL(allLvalues(idx)) = cmap(idx, :);
end

% Distinct marker SHAPE per L too, not just color -- so if two lines sit
% on top of each other (e.g. two near-zero-BER curves), you can still
% tell them apart instead of one fully hiding the other.
markerShapes = {'o', 's', '^', 'd', 'v', 'p', 'h'};
markerForL = containers.Map('KeyType', 'double', 'ValueType', 'any');
for idx = 1 : numel(allLvalues)
    markerForL(allLvalues(idx)) = markerShapes{mod(idx-1, numel(markerShapes)) + 1};
end

%% BER vs SNR, all effective L values on one axes, one figure per decoder
for dec = 1 : numDecoders
    figure('Name', sprintf('BER vs SNR by effective L (%s)', decoderNames{dec}), ...
        'Color', 'w', 'Position', [100 100 800 500]);
    hold on;

    % Plotted in REVERSE (largest to smallest) so the smallest L -- most
    % likely to sit at BER=0 -- is drawn LAST and ends up on top/visible
    % instead of buried under the others.
    for Lval = fliplr(Lvalues_toPlot)
        rows = (Lall == Lval);
        if ~any(rows), continue; end
        [snrSorted, order] = sort(snrDbAll(rows));
        berRows = berAll(rows, dec);
        berSorted = berRows(order);
        plot(snrSorted, berSorted * 100, '-', 'Marker', markerForL(Lval), ...
            'LineWidth', 2, 'MarkerSize', 7, ...
            'Color', colorForL(Lval), 'MarkerFaceColor', colorForL(Lval), ...
            'DisplayName', sprintf('L = %d', Lval));
    end
    grid on; box on;
    xlabel('Relative SNR [dB]  (10 log_{10}(power scale))');
    ylabel('BER [%]');
    title(sprintf('%s: BER vs SNR, by effective guard length', decoderNames{dec}));
    legend('Location', 'northeastoutside');
    set(gca, 'FontSize', 11);

    saveas(gcf, sprintf('BER_vs_SNR_by_effectiveL_%s.png', decoderNames{dec}));
end