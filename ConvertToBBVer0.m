function[ SignalBB, FiltMem, NextPhase, FiltMemDec] = ConvertToBBVer0( ...
                                           Signal, Fcarrier, Fs, Factor,...
                                           FiltCoeff, FiltMem, Phase)               
% function[ SignalBB, FiltMem, NextPhase, FiltMemDec ] = ConvertToBBVer0( ...
%                                            Signal, Fcarrier, Fs, Factor,...
%                                            FiltCoeff, FiltMem, Phase)
% 
% Description:
%           This function convert a crossband signal to a baseband signal. 
%
% inputs:
%           Signal      - The signal to be converted.
%           Fcarrier    - carrier frequency of the sent signal.
%           Fs          - sampling frequency.
%           Factor      - Decimation Factor.
%           FiltMem     - Filter memory vector of length(FiltCoeff)-1...  
%                         if zero - there is no memory from the last time.
%           Phase       - carrier phase from the last call of the function.
%           FiltCoeff   - The LPF coefficients.
%                             
% Output:                         
%           SignalBB    - The signal in base band. In the first call of the
%                         function, it has a delay of length(FiltCoeff)/2 samples 
%           FiltMem     - Filter memory for the next call.
%           NextPhase   - Carrier phase to the next call of the function.
%           FiltMemDec  - Filter memory after decimation - Concatenate it to ...
%                                'SignalBB' in the last call of the function.
% 
% Comments: Modified by Kobi Bucris 10.7.07
%           Updated by Shlomi Museri 15.9.11
% 
% 

% Set default values.
if nargin < 6
    FiltMem = zeros(1,length(FiltCoeff)-1);
end
if nargin < 7
    Phase = 0;
end

DataLen =  length(Signal); % The signal length.
% Warning message.
if mod(DataLen, Factor)
    warning('Data length is not a multiple of the decimation factor,there could be some inaccuracies in the adges');
end

time = 0:1/Fs:(DataLen-1)/Fs; % Time vector.
% Calculate the NextPhase.
NextPhase = (time(end) + 1/Fs ) * 2 * pi * Fcarrier + Phase; 
NextPhase  = mod(NextPhase, 2*pi);

SignalBB = Signal .* exp(-1i*(2 * pi * Fcarrier * time + Phase)) ; % Shift the signal to BaseBand.
% Filtering and resampling the baseband signal.
[SignalBB,FiltMem] = filter( FiltCoeff, 1, SignalBB ,FiltMem);
SignalBB = SignalBB(1:Factor:end);  % Decimation.
FiltMemDec = FiltMem(1:Factor:end); % Decimation.
FiltMem = FiltMem.';       % Convert to row vector.
FiltMemDec = FiltMemDec.'; % Convert to row vector.