% LoRaDetector_FFTsegments.m
% ---------------------------------------------
% This script generates a LoRa signal, passes it through a wireless channel,
% and visualizes the FFT output of each segment of length M (no preamble used).
% ---------------------------------------------

close all; clear all; clc;

plotDuringRun = false; % Set to true to enable plotting during SNR sweep

% SNR sweep parameters
SNRsweep = -20:1:-10; % dB

minSNRs = 4.4;%4:0.2:6; % Set your desired minSNR threshold(s) here
testVals = minSNRs;
misdetectRate = zeros(length(testVals), length(SNRsweep));
falseAlarmRate = zeros(length(testVals), length(SNRsweep));
numTrials = 1000; 

% Parameters
SF = 7;           % Spreading Factor
BW = 125e3;       % Bandwidth (Hz)
Fs = BW;          % Original sampling frequency
M = 2^SF;         % Number of frequency bins (segment length)
Nsymbols = 100;   % Number of symbols to transmit
Nbits = Nsymbols * SF;

% Channel params
timingOffset = M/16 + 0.5; % Test at this specific timing offset (fractional sample)
maxDelay = 0;%timingOffset;    % For clarity, use this as the delay in the channel
max_Fo = BW;                % No frequency offset
dF = BW/M;

for testIdx = 1:length(testVals)
    for snrIdx = 1:length(SNRsweep)
        SNR = SNRsweep(snrIdx);
        misdetects = zeros(1, numTrials);
        falseAlarms = zeros(1, numTrials);
        for trial = 1:numTrials

            % Generate random data bits
            TxDataBits = randi([0 1], Nbits, 1);

            % Modulate data (no preamble)
            modulatedSignal = CSSmod(reshape(TxDataBits, SF, []).', SF, BW, 1, true, false);
            modulatedSignal = modulatedSignal(:); % Ensure column vector


            % Apply timing offset (fractional delay) before channel
            delaySamples = timingOffset;
            % % Use linear interpolation for fractional delay
            t = (0:length(modulatedSignal)-1)';
            t_delayed = t - delaySamples;
            modulatedSignalDelayed = interp1(t, modulatedSignal, t_delayed, 'spline', 0);

            % Pass through wireless channel (no additional delay)
            [RxSignal, ~, ~] = wireless_channel(modulatedSignalDelayed, Fs, 0, max_Fo, SNR, dF);

            detected = lora_detect_signal(RxSignal, SF, BW, testVals(testIdx), M, plotDuringRun);
            misdetects(trial) = ~detected;

            % Noise-only test for false alarm
            noiseOnly = randn(length(RxSignal),1) + 1j*randn(length(RxSignal),1);

            detectedNoise = lora_detect_signal(noiseOnly, SF, BW, testVals(testIdx), M, plotDuringRun);
            falseAlarms(trial) = detectedNoise;
        end
        misdetectRate(testIdx, snrIdx) = mean(misdetects);
        falseAlarmRate(testIdx, snrIdx) = mean(falseAlarms);
        fprintf('TestVal = %.2f, SNR = %d dB, Misdetect Rate = %.2f, False Alarm = %.2f\n', testVals(testIdx), SNR, misdetectRate(testIdx, snrIdx), falseAlarmRate(testIdx, snrIdx));
    end
end

% Plot FAMD vs SNR

figure;
hold on;
colors = lines(length(testVals));
for testIdx = 1:length(testVals)
    figure;
    plot(SNRsweep, misdetectRate(testIdx,:), '-o', 'Color', colors(testIdx,:), 'LineWidth', 1.5, 'DisplayName', ['MD Thresh=' num2str(testVals(testIdx))]);
    hold on;
    plot(SNRsweep, falseAlarmRate(testIdx,:), '-x', 'Color', colors(testIdx,:), 'LineWidth', 1.5, 'DisplayName', ['FA Thresh=' num2str(testVals(testIdx))]);
    xlabel('SNR (dB)'); ylabel('Misdetect and False Alarm Rate');
    title('Misdetect Rate vs SNR for Different Thresholds');
    legend('show');
    grid on;
end



% --- Helper function for LoRa signal detection over all segments ---

function detected = lora_detect_signal(RxSignal, SF, BW, threshold, M, plotDuringRun)
    numSegments = floor(length(RxSignal)/M);
    fft_centered_all = zeros(M, numSegments);
    bitWord = de2bi(0, SF, 'left-msb');
    downchirp = CSSmod(bitWord, SF, BW, -1, true, false);

    for segIdx = 1:numSegments
        segment = RxSignal((segIdx-1)*M+1 : segIdx*M);
        dechirped = segment .* downchirp;
        fft_raw = fftshift(fft(dechirped));
        [~, peakIdx] = max(abs(fft_raw));
        shift_amt = (M/2 + 1) - peakIdx;
        fft_centered = circshift(fft_raw, shift_amt);
        fft_centered_all(:, segIdx) = fft_centered;
        if plotDuringRun
            figure(1001); hold on;
            plot(abs(fft_centered), 'Color', [0.7 0.7 0.7 0.5]);
        end
    end

    fft_combined = mean(abs(fft_centered_all),2);

    if plotDuringRun
        figure(1001); hold on;
        plot(fft_combined, 'k', 'LineWidth', 2, 'DisplayName', 'Mean Combined');
        title('All Centered FFTs (Overlayed) and Mean');
        xlabel('Frequency Bin (centered)'); ylabel('Magnitude');
        grid on;
    end

    peakVal = fft_combined(M/2+1);
    wrapIdx = mod((M/2-2:M/2+4)-1, M) + 1;
    noiseVals = fft_combined;
    noiseVals(wrapIdx) = [];
    noiseMean = mean(noiseVals);
    metric = 10*log10(peakVal / noiseMean);
    detected = metric > threshold;
    if plotDuringRun
        fprintf('Combined FFT: Peak=%.2f, NoiseMean=%.2f, Metric=%.2f, Threshold=%.2f, Detected=%d\n', peakVal, noiseMean, metric, threshold, detected);
        title(['PPSNR = ' num2str(metric) '[dB]']);
    end
end

% --- Helper function for segment metric calculation ---
% (Optional) Segment-level metric calculation function
% Not used in main detection, but useful for debugging or alternate logic
function [metric, detectedSeg, fft_out, peakIdx] = lora_segment_detection(segment, SF, BW, threshold)
    bitWord = de2bi(0, SF, 'left-msb');
    downchirp = CSSmod(bitWord, SF, BW, -1, true, false);
    dechirped = segment .* downchirp;
    fft_out = abs(fftshift(fft(dechirped)));
    [peakVal, peakIdx] = max(fft_out);
    Nfft = length(fft_out);
    wrapIdx = mod((peakIdx-4:peakIdx+3)-1, Nfft) + 1;
    noiseVals = fft_out;
    noiseVals(wrapIdx) = [];
    noiseMean = mean(noiseVals);
    metric = 10*log10(peakVal / noiseMean);
    detectedSeg = metric > threshold;
end


