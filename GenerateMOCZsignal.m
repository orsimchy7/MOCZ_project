clear all;
close all;

% -- Params --
fileName = 'MOCZ_tx.mat';
B = 200; %num of sequences in a message. =20
L = 5; % additional taps that the channel conv adds
K = 12 ; %X(z) polynomial order. there are K+1 coefficents
Tsym = 1e-3;
fc = 15e3;
Fs = 200e3;
lambda = 1;
figFlag = 0;
betha = 1;
BW = 1 / Tsym; 
BWn = (1+betha)*BW;

% adding zadOffChu for sync:
N = 63;     % Sequence length
u = 2;     % Root index (gcd(u,N)=1)
zadoffChuPair = [u, N];


% P is a binary packet. there are B messages (columns) and K bits in each
% message
P = randi([0,1], K, B); 
signal_pb_total = [];

% -- Huffman BMOCZ --
Kidxs = 1:K;
R = sqrt(1+2*lambda*sin(pi/K));
theta_c = ((2*pi) * (Kidxs / K))';

simParams = struct('B', B, 'K', K, 'Tsym', Tsym, 'fc', fc, 'Fs', Fs, ...
    'lambda', lambda, 'L', L, 'figFlag', figFlag, 'betha', betha, ...
    'zadoffChuPair', zadoffChuPair, 'R', R, 'theta_c', theta_c);

% for sync:
zc = zadoffChuSeq(u, N);
% Initialize an empty vector for all baseband symbols
symbols_total = [];
% Add Zadoff-Chu sequence and its guard to the beginning
zc_with_guard = [zc; zeros(L, 1)];

symbols_total = [symbols_total; zc_with_guard];

% Process each message in the packet
for b = 1:B
    % Get current message bits
    M = P(:, b);
    
    % Generate Alphas (Roots)
    alphas = ((1-M)*(R^(-1)) + M * R) .* exp(1i*theta_c);
    
    % Generate Coefficients
    x_un = flip(poly(alphas))'; %finding K+1 coefficients and flip - now the free coefficent is first
    x = (x_un / norm(x_un)) * sqrt(K+1); % normalize to sqrt(K+1). now total energy of signal will be K+1
    
    x_with_guard = [x; zeros(L, 1)]; %Guard of zeros (channel of L taps)
    
    symbols_total = [symbols_total; x_with_guard];

end

[signal_pb_total, ~] = pulseSHP(symbols_total, simParams, 'modulate');

% Save to .mat file
fileName = sprintf('L%dB%dK%dN%du%d', L, B, K, N, u);
% path = fullfile("C:\Users\SHOHAMM\Desktop\projB\source_code\MOCZ_pool_updated_17.2\signalsPool", fileName);
%path = fullfile("C:\Users\user\OneDrive - Technion\Desktop\פרויקט ב 5.3.26", fileName);

%save in the script's location:
scriptPath = fileparts(mfilename('fullpath'));
path = fullfile(scriptPath, fileName);

save(path, 'signal_pb_total', 'P', 'simParams');
fprintf('MOCZ signal with B=%d messages saved to %s\n', B, fileName);