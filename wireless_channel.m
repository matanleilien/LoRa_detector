function [y,delay,Foffset] = wireless_channel(x, fs, max_delay, max_freq_offset, snr_dB,dF)
  % x: input signal asusmed to be column
  % fs: sampling frequency (Hz)
  % max_delay: maximum delay in samples
  % max_freq_offset: maximum frequency offset in Hz
  % snr_dB: signal-to-noise ratio in dB

  if size(x,1) == 1 %make x column if it's not
    x = x.';
  end

  % Random delay
  delay = randi([0, max_delay]);
  x_delayed = [zeros(delay,1);x;zeros(floor(length(x)/4),1)];
  N = length(x_delayed);

  % Random frequency offset
  Foffset = (2*rand - 1) * max_freq_offset;  % Uniform in [-max, max]

  FoInt = floor(Foffset / dF) * dF;
  FoFrac = Foffset - FoInt;

  % Enable fractional freq offset
  Foffset = FoInt + FoFrac;
  % disp(['freq offset is : ' num2str(Foffset/dF) ' times dF, frac offset = ' num2str(FoFrac) ' Hz']);

  t = (0:N-1)' / fs;
  freq_shift = exp(1j * 2 * pi * Foffset * t);
  x_freq_shifted = x_delayed .* freq_shift;

  % Add AWGN
  signal_power = mean(abs(x_freq_shifted).^2);
  noise_power = signal_power / (10^(snr_dB/10));
  noise = sqrt(noise_power/2) * (randn(size(x_freq_shifted)) + 1j*randn(size(x_freq_shifted)));
  y = x_freq_shifted + noise;
end

