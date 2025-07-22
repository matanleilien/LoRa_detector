function cssSamples = CSSmod(dataBits, SF, BW, chirpSign, doFM,verbose)
% CSSmod - Chirp Spread Spectrum (CSS) modulation
%
%   cssSamples = CSSmod(dataBits, SF, BW, chirpSign, verbose, doFM)
%
%   Inputs:
%     dataBits  : Vector of input bits to modulate
%     SF        : Spreading Factor (number of bits per symbol)
%     BW        : Bandwidth (Hz)
%     chirpSign : 1 for up-chirp, -1 for down-chirp (default is 1)
%     verbose   : (optional) true to plot, false to suppress plots (default: false)
%     doFM      : (optional) true to output FM modulated signal (default: true)
%
%   Output:
%     cssSamples: Vector of complex CSS modulated samples (if doFM true)
%                or matrix of instantaneous frequency values (if doFM false)
%
%   The function pads dataBits if needed, splits into symbols, and generates
%   a CSS modulated signal for each symbol. It can also plot the instantaneous
%   frequency and the spectrogram for each symbol if verbose is true.

if nargin < 4
    chirpSign = 1;
end
if nargin < 5
    doFM = true;
end
if nargin < 6
    verbose = false;
end

Fs = BW; % Sampling frequency, set to bandwidth for simplicity
bitsPerSymbol = SF; % Number of bits per symbol
M = 2^SF;
sps = M; % Samples per symbol

dataBits = dataBits(:)'; % Ensure row vector for bit grouping

% Pad dataBits with zeros if not a whole multiple of bitsPerSymbol
numBits = length(dataBits);
remainder = mod(numBits, bitsPerSymbol);
if remainder ~= 0
    padLength = bitsPerSymbol - remainder;
    dataBits = [dataBits, zeros(1, padLength)];
end
Nsymbols = length(dataBits) / bitsPerSymbol; % Number of symbols

%     k = (0:M-1);              % k vector
t = (0:1/Fs:(M - 1)/Fs); % Time vector for one symbol

df = BW / M; % Frequency step

cssSamples = [];


for curSymbol = 1:Nsymbols
    bits = dataBits((curSymbol-1)*bitsPerSymbol+1 : curSymbol*bitsPerSymbol);
    % Convert bits to integer symbol (MSB first)
    m = bi2de(bits, 'left-msb');
    fm = -BW/2 + m*df;
    fi =fm + chirpSign*t*BW*df;


    if verbose
        % Plot instantaneous frequency
        figure;
        plot(t, fi);
        title(["fi for symbol ", num2str(curSymbol), " (chirpSign = ", num2str(chirpSign), ")"]);
        xlabel("Time (s)"); ylabel("fi");
    end
    if doFM
        % FM Modulation
        phase = 2*pi*cumtrapz(t,fi);

        fmSignal = (exp(1i*phase)).'; % Complex baseband CSS signal
        cssSamples = [cssSamples; fmSignal]; % Concatenate as column vector
        if verbose
            plotCSSSpectrogram(fmSignal, Fs, M, curSymbol,-20);
        end
    else
        cssSamples = [cssSamples; fi(:)];
    end
end

end
