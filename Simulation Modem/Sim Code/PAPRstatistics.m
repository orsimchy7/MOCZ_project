% Validation script - CCDF of PAPR for K=5 and K=12 side-by-side
clear;

% Define the K values to evaluate
K_vals = [5, 12];

% Create a wider figure to accommodate side-by-side subplots
figure('Position', [100, 100, 1000, 450]);

for idx = 1:length(K_vals)
    K = K_vals(idx);
    lambda = 1;
    Kidxs = 0:K-1;
    R = sqrt(1+2*lambda*sin(pi/K));
    theta_c = ((2*pi) * (Kidxs / K))';
    
    bits = ff2n(K);
    num_combinations = 2^K;
    PAPR = zeros(num_combinations, 1);
    
    for ii = 1 : num_combinations
        % Calculate the PAPR for the current bit sequence:
        currentBits = bits(ii, :);
        M = currentBits';
        
        % 1. Generate Alphas (Roots)
        alphas = ((1-M)*(R^(-1)) + M * R) .* exp(1i*theta_c);
        
        % 2. Generate Coefficients
        x_un = flip(poly(alphas)).'; % finding K+1 coefficients and flip
        x = (x_un / norm(x_un)) * sqrt(K+1); % normalize to sqrt(K+1)
        
        % Calculate PAPR in dB
        PAPR(ii) = 10 * log10(max(abs(x).^2) / mean(abs(x).^2));
    end
    
    % Calculate probability P(PAPR > threshold)
    papr_sorted = sort(PAPR); 
    N = length(PAPR);
    prob = (N:-1:1) / N; 
    
    % Plot on the corresponding subplot
    subplot(1, 2, idx);
    semilogy(papr_sorted, prob, 'LineWidth', 2);
    grid on;
    xlabel('PAPR_0 [dB]');
    ylabel('Prob(PAPR > PAPR_0)');
    title(sprintf('CCDF of PAPR (K = %d)', K));
    
    % Optional: Set consistent X-limits to make the visual shift obvious
    % xlim([2 12]); 
end

%%
% Validation script - PAPR at specific probability threshold vs K

% --- Configuration ---
target_prob = 0.01; % Set to 0.01 for 1%, 0.10 for 10%, 0.001 for 0.1%, etc.
K_vals = 5:12;      % Block sizes to evaluate
% ---------------------

papr_at_target = zeros(size(K_vals));

for idx = 1:length(K_vals)
    K = K_vals(idx);
    lambda = 1;
    Kidxs = 0:K-1;
    R = sqrt(1+2*lambda*sin(pi/K));
    theta_c = ((2*pi) * (Kidxs / K))';
    
    num_combinations = 2^K;
    bits = ff2n(K);
    PAPR = zeros(num_combinations, 1);
    
    for ii = 1 : num_combinations
        M = bits(ii, :)';
        alphas = ((1-M)*(R^(-1)) + M * R) .* exp(1i*theta_c);
        x_un = flip(poly(alphas)).';
        x = (x_un / norm(x_un)) * sqrt(K+1);
        PAPR(ii) = 10 * log10(max(abs(x).^2) / mean(abs(x).^2));
    end
    
    % Sort and calculate probability
    papr_sorted = sort(PAPR);
    prob = (num_combinations:-1:1) / num_combinations;
    
    % Use linear interpolation to extract the exact PAPR at the target probability.
    % This eliminates the "jaggedness" caused by jumping between discrete array indices.
    papr_at_target(idx) = interp1(prob, papr_sorted, target_prob, 'linear');
end

% Plotting the relationship
figure(3);
plot(K_vals, papr_at_target, '-o', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('Number of Bits (K)');
ylabel(sprintf('PAPR at %.1f%% Probability [dB]', target_prob * 100));
title(sprintf('%.1f%% CCDF PAPR vs. Block Size (K)', target_prob * 100));