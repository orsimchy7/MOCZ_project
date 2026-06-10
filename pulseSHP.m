function [signal_pb, x_decoded, group_delay] = pulseSHP(x, simParams, mode)
%Description:
% The function had two operational modes: "modulate" / "demodulate"
%   modulate: the function gets coefficients to transmit and apply
%   raised-cosine pulse shaping to create BB signal. later it creates and
%   returns the passBand signal
    
    % -- Params --
    L = simParams.L;
    K = simParams.K;
    N = K + L;
    Tsym = simParams.Tsym;
    fc = simParams.fc;
    Fs = simParams.Fs;
    beta = simParams.betha;
 
    sps = Tsym * Fs; %samples per symbol
    span = 6;
    h_coeff = rcosdesign(beta, span, sps); %taps num is span*sps + 1
    normalization_factor = 1;
    group_delay = span * sps;
    t = (0:1/Fs:Tsym*(N+1 +group_delay/sps) - 1/Fs)';


    if strcmp(mode, 'modulate')

        x_upsamp = upsample(x, sps);
        x_upsamp = [x_upsamp; zeros(group_delay,1)];
        signal_bb = filter(h_coeff, normalization_factor, x_upsamp);
        signal_bb= x_upsamp;
        t = (0:length(signal_bb)-1)' / Fs;
        x_decoded = [];
        signal_pb = sqrt(2) * real(signal_bb.*exp(1i*2*pi*fc*t));
        
    elseif strcmp(mode, 'demodulate')
        % -- Decoding raised cosing --
        t = (0:1/Fs:Tsym*(K + 1 + (group_delay/sps)) - 1/Fs)';
        t = (0:length(x)-1)' / Fs;
        mixed_x = x.*exp(-1i*2*pi*fc*t);
        x_bb_filtered = filter(h_coeff, normalization_factor, mixed_x);
        
        % Downsampling
        % We start sampling after the delay and then take 1 sample every 'sps'
        x_decoded = sqrt(2) * x_bb_filtered(group_delay + 1 : sps : end);
        signal_pb = [];
    end
end