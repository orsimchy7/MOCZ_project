function [bits_rx_greedy] = Greedy(x_vec_decoded,simParams)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    K = simParams.K;
    R = simParams.R;
    N = simParams.N; %L+K
    V_bank = simParams.V_bank;
    
    noisy_signal = (x_vec_decoded(1:end-1))'; % take only L-1 guard out of L!
    y = noisy_signal(:); % Ensure y is a column vector of size Nx1. no need for last guard!

    %Initialization with DiZeT
    bits_rx_greedy = zeros(K,1);
    zero_1 = R;
    conj_zero_1 = (1/R);
    Y_zk = abs(polyval(fliplr(noisy_signal), zero_1));
    Y_zk_conj = abs(polyval(fliplr(noisy_signal), conj_zero_1));
    if Y_zk > (R^(N))*Y_zk_conj
        bits_rx_greedy(1) = 0;
        V_temp = V_bank(:,K+1);
    else
        bits_rx_greedy(1) = 1;
        V_temp = V_bank(:,1);
    end    
    mask = zeros(1,2*K);
    mask(1) = 1;
    mask(K+1) = 1;
    P_i = V_temp*inv(V_temp'*V_temp)*V_temp'; %NxN
    P_i_orth = eye(N,N)-P_i;
    
    % greedy selection loop
    for t = 2:K
        min_cost = inf;
        min_l = 0;
        for l = 1:2*K
            if mask(l) == 1
                continue
            end
            u_i = (P_i_orth*V_bank(:,l))/norm((P_i_orth*V_bank(:,l)));
            projection = (abs(u_i'*y))^2;
            if projection < min_cost
                min_cost = projection;
                min_l = l;
            end
        end
        V_temp = [V_temp V_bank(:,min_l)];
        if min_l > K
            k_idx = min_l - K;
            mask(min_l) = 1;
            mask(min_l - K) = 1;
            bit_val = 0;
        else
            k_idx = min_l;
            mask(min_l) = 1;
            mask(min_l+K) = 1;
            bit_val = 1;
        end
        bits_rx_greedy(k_idx) = bit_val;
    
        P_i = V_temp*inv((V_temp)'*V_temp)*V_temp';    
        P_i_orth = eye(N,N)-P_i;
    end
end