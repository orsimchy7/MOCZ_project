function [bits_rx_ML, W_mat] = ML(x_vec_decoded,simParams)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


    K = simParams.K;
    noisy_signal = (x_vec_decoded(1:end-1)).'; % take only L-1 guard out of L!
    y = noisy_signal(:); % Ensure y is a column vector of size Nx1. no need for last guard!


    combinations = simParams.combinations;
    Cinv_all = simParams.Cinv_all;
    V_all = simParams.V_all;

    % decoding with ML
    min_cost = inf;
    bits_rx_ML = zeros(K, 1);
    
    for idx = 1:2^K
        V = V_all(:,:,idx);          % Vandermonde matrix 
        Cinv = Cinv_all(:,:,idx);    % Inverse sqrt Gram matrix (K×K)

        Vy = V' * y;                  % (K×1)
        projection = Cinv * Vy;      % (K×1)
        cost = norm(projection)^2;

        if cost < min_cost
            min_cost = cost;
            bits_rx_ML = combinations(idx, :).';
            W_mat = Cinv; %diagonal or not?
        end
    end      


end