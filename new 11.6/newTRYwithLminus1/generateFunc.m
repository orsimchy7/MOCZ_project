function [P, signal_pb_total, x_bb, group_delay, simParams] = generateFunc(simParams)
% -- Params --
% MOCZ
L = simParams.L;
K = simParams.K;
% DEBUG:
% K=2;

% Check for image coding flag (defaults to 0 if not provided)
if isfield(simParams, 'image_coding')
    image_coding = simParams.image_coding;
else
    image_coding = 0;
end

image_power_test = 0; %if = 1 -> will decrease power of signal
image_snr_test = 0; %if = 1 -> will add awgn

% zadOffChu (sync):
N = simParams.zadoffChuPair(2);
u = simParams.zadoffChuPair(1);

%% Packet Configuration (Image vs. Random)
if image_coding == 1
    % --- IMAGE CODING MODE ---
    
    % Determine which image file to read (checking common names)
    target_file = '';
    if isfield(simParams, 'image_file') && isfile(simParams.image_file)
        target_file = simParams.image_file;
    elseif isfile('my_image.jpg')
        target_file = 'my_image.jpg';
    end
    
    if ~isempty(target_file)
        % Read the original image
        img_raw = imread(target_file);
        
        % Convert logical matrix to double to avoid downstream function crashes
        if islogical(img_raw)
            img_raw = double(img_raw);
        end
        
        % Resize to exactly 256 x 256 pixels
        img_resized = imresize(img_raw, [256, 256]);
        
        % Convert to grayscale if it is a color image (RGB)
        if size(img_resized, 3) == 3
            img_gray = rgb2gray(img_resized);
        else
            img_gray = img_resized;
        end
        
        % Convert to binary (pure black and white - logical 0 and 1)
        img_binary = imbinarize(img_gray);
    else
        % Fallback: Generate a default 256x256 geometric cross pattern if no file exists
        img_binary = zeros(256, 256);
        img_binary(88:168, 88:168) = 1;  % Center square
        img_binary(32:224, 116:140) = 1; % Vertical bar
        img_binary(116:140, 32:224) = 1; % Horizontal bar
        fprintf('Warning: Image file not found. Using default 256x256 cross pattern.\n');
    end
    
    % Serialize the 256x256 image into a 1D bitstream (65,536 bits)
    bit_stream = double(img_binary(:));
    
    % Calculate and update B dynamically to fit the image exactly
    B = length(bit_stream) / K; 
    simParams.B = B; 
    
    % Shape into the K x B packet matrix
    P = reshape(bit_stream, [K, B]);
    simParams.P = P;
    simParams.power = 100;
    
else
    % --- RANDOM PERMUTATIONS MODE  ---
    B = simParams.B;
    P = zeros(K, B);
    
    usedIdx = simParams.usedIdx;
    bits = ff2n(K);
    for jj = 1 : B
        idx = randi([1, 2^K]);
        while ~all(usedIdx) && usedIdx(idx)
            idx = randi([1, 2^K]); 
        end
        usedIdx(idx) = true; 
        P(:, jj) = bits(idx, :)'; 
    end
    simParams.usedIdx = usedIdx;
end

% P = ones(K,B); %DEBUG!!!
% P = zeros(K,B); %DEBUG!!!
simParams.P = P; %we save the packet's bits. so we can check error rate of data.

signal_pb_total = [];

%% -- Huffman BMOCZ --
R = simParams.R;
theta_c = simParams.theta_c;

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
    x_un = flip(poly(alphas)).'; %finding K+1 coefficients and flip - now the free coefficent is first
    x = (x_un / norm(x_un)) * sqrt(K+1); % normalize to sqrt(K+1). now total energy of signal will be K+1
    
    x_with_guard = [x; zeros(L-1, 1)]; %Guard of zeros (channel of L taps)
    
    symbols_total = [symbols_total; x_with_guard];
end


[signal_pb_total, ~, group_delay] = pulseSHP(symbols_total, simParams, 'modulate');
% no pulse shaping:
%signal_pb_total = symbols_total;
% group_delay = 1;
M_sig = max(signal_pb_total, [], "all");
m_sig = min(signal_pb_total, [], "all");
a = simParams.a;
b = simParams.b;
stretched_sig = a + (signal_pb_total - m_sig) * (b - a)/(M_sig - m_sig);
% stretched_sig = (signal_pb_total/M_sig)*b;
x_bb = stretched_sig;
signal_pb_total = stretched_sig;

simParams.P = P;

%Change image power
if image_coding == 1 && image_power_test ==1
    power = 0.2; %can be in range [0 1]
    amp_scale = sqrt(power);
    x_bb = x_bb * amp_scale;
    signal_pb_total = signal_pb_total * amp_scale;
    simParams.power = power;
end

if image_coding == 1 && image_snr_test ==1
    SNR = -8; %[dB]
    x_bb = awgn(x_bb, SNR, 'measured');
    signal_pb_total = awgn(signal_pb_total, SNR, 'measured');
    simParams.testSNR = SNR;
end



end