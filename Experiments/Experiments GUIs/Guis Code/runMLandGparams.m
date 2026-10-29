function [simParamsNew] = runMLandGparams(simParams, K, L)
    simParamsNew = simParams;
    R = simParams.R;
    theta_c = simParams.theta_c;
    % 3. Greedy decoder
    N = K + L ;   % Convolution length
    zeros_bank = zeros(1,2*K); %2K optional zeros of polynom
    zeros_bank(1,1:K) = R*exp(1i * theta_c(1:K)); %ONES
    zeros_bank(1,(K+1):2*K) = (1/R)*exp(1i * theta_c(1:K)); %ZEROS
    V_bank = zeros(N, 2*K); %all Vj vectors
    for col = 1:2*K
        V_bank(:, col) = conj(zeros_bank(col)).^(0:N-1);
    end
    
    % 4. ML decoder
    % Precompute all 2^K binary combinations and their corresponding zero vectors
    combinations = dec2bin(0:2^K-1) - '0';   %matrix of all 2^k possible messages
    V_all = zeros(N,K,2^K); %all possible vandermonde matrixes
    
    alpha_candidates = zeros(2^K, K); %all possible alphas - mapping
    for idx = 1:2^K
        for k = 1:K
            if combinations(idx, k) == 1
                alpha_candidates(idx, k) = R * exp(1i * theta_c(k));
            else
                alpha_candidates(idx, k) = (1/R) * exp(1i * theta_c(k));
            end
        end
        
        alpha = alpha_candidates(idx, :);
        V = zeros(N, K);
        for col = 1:K
            V(:, col) = conj(alpha(col)).^(0:N-1);
        end
    
        V_all(:,:,idx) = V;
    
        C = V' * V;
        [U, S_diag, ~] = svd(C);
        Cinv_all(:,:,idx) = U * diag(1 ./ sqrt(diag(S_diag))) * U';
        % Cinv_all(:,:,idx) = inv(sqrtm(C));
    end

    simParamsNew.N = N;
    simParamsNew.V_bank = V_bank;
    simParamsNew.combinations = combinations;
    simParamsNew.Cinv_all = Cinv_all;
    simParamsNew.V_all = V_all;
end