folderPath = 'C:\Users\Talkid\Desktop\MOCZ_Pool\Recordings_pool_7_6_26\fiveVmax\MOCZ_packets\L_changes';
files = dir(folderPath);
files = files(3:end -1, :);
M = length(files) / 5 - 2; % assumes try = 5;
L = zeros(length(files) / 5 ,1);
berL = zeros(length(files) / 5 ,1);
for k = 1 : length(files) / 5
   for j = 1 : 5
        idx = 5 * (k-1) + j;
        basefileName = files(idx).name;
        fullfileName = fullfile(folderPath, basefileName);
        fprintf('now: %s\n', basefileName);
        load(fullfileName);
        L(k) = currentParams.L;
        DataMatrix = [RecordedData.Variables];
        signalData1 = DataMatrix(:, 1);
        signalData2 = DataMatrix(:, 2);
        [~, error_rate1, ~] = decode_mocz(signalData1, currentParams);
        [~, error_rate2, ~] = decode_mocz(signalData2, currentParams);
        errorR = (error_rate1 + error_rate2) * 0.5;
        berL(k) = berL(k) + 0.2 * errorR;
   end

end