function [message_rec] = DiZeT(R, theta_c, y, K)
    %UNTITLED2 Summary of this function goes here
    %   x_rec are reconstructed polinom coefficients of x
    % UPDATE: convert y so dizet will be symetric to ML (decode from the
    % same symbols number) + treat edge case of conv with 0 coeff:
    if length(y) > K + 1
        y = y(1:end-1);
    end
    Lt = length(y);
    message_rec = zeros(K,1);
    
    for t = 1:length(theta_c)
        x_rec_flip = flip(y.');
        inVal = polyval(x_rec_flip, R^(-1) * exp(1i*theta_c(t)));
        outVal = polyval(x_rec_flip, R * exp(1i*theta_c(t)));
        if abs(outVal) < R^(Lt-1)*abs(inVal)
            message_rec(t) = 1;
        end
    end
end