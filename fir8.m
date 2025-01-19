% Frequency domain analysis of the original signal
nfft = length(audio); % FFT length
f = (0:nfft-1)*(ws/nfft); % Frequency vector
originalFFT = abs(fft(audio, nfft)); % Magnitude of FFT

% Plot the original signal frequency spectrum
figure;
plot(f, originalFFT);
title('Original Signal Frequency Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% Normalize the input signal
input_signal = input_signal / max(abs(input_signal));

% Define cutoff frequencies
low_cutoff = 300;    % Lower cutoff frequency in Hz
high_cutoff = 1500;  % Upper cutoff frequency in Hz

% Check if cutoff frequencies are valid
nyquist_frequency = fs / 2; % Nyquist frequency
if high_cutoff >= nyquist_frequency
    error('High cutoff frequency must be less than Nyquist frequency. Please adjust.');
end

% Normalize the cutoff frequencies
normalized_cutoff = [low_cutoff, high_cutoff] / nyquist_frequency;

% Design FIR bandpass filter with Hamming window for N=8
N = 8; % Filter order
b = fir1(N, normalized_cutoff, 'bandpass', hamming(N+1));

% Apply the filter to the signal
filtered_signal = filtfilt(b, 1, input_signal); % Zero-phase filtering

% Save the filtered audio file
audiowrite('filtered_audio_FIR_Hamming_N8_300to1500.wav', filtered_signal, fs);

% Plot (d): Frequency response of the filter |H(ω)|, |H(ω)|_dB, ∠H(ω)
[H, W] = freqz(b, 1, 1024, fs);
figure;
subplot(3, 1, 1);
plot(W, abs(H), 'b', 'LineWidth', 1.5);
title('Frequency Response |H(ω)|');
xlabel('Frequency (Hz)');
ylabel('|H(ω)|');
grid on;

subplot(3, 1, 2);
plot(W, 20*log10(abs(H)), 'r', 'LineWidth', 1.5);
title('Frequency Response in dB |H(ω)|_dB');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;

subplot(3, 1, 3);
plot(W, angle(H) * 180 / pi, 'g', 'LineWidth', 1.5);
title('Phase Response ∠H(ω)');
xlabel('Frequency (Hz)');
ylabel('Phase (degrees)');
grid on;

% Plot (e): Pole-zero plot of the filter H(z)
figure;
zplane(b, 1);
title('Pole-Zero Plot of FIR Filter H(z)');

% Plot (f): Impulse response h[n]
figure;
stem(b, 'filled');
title('Filter Impulse Response h[n]');
xlabel('Sample Index (n)');
ylabel('Amplitude');

% Plot (g): Frequency components after filtering |Y(ω)|
output_fft = fft(filtered_signal);
figure;
plot(frequencies(1:floor(length(frequencies)/2)), abs(output_fft(1:floor(length(frequencies)/2))));
title('Frequency Spectrum After Filtering |Y(ω)|');
xlabel('Frequency (Hz)');
ylabel('|Y(ω)|');

% Plot (h): Output signal in time domain y[n]
figure;
plot(filtered_signal);
title('Output Signal in Time Domain (y[n])');
xlabel('Sample Index');
ylabel('Amplitude');

% Play the filtered audio
disp('Playing filtered audio...');
sound(filtered_signal, fs);

% Display a summary
disp('Filtered audio saved as: filtered_audio_FIR_Hamming_N8_300to1500.wav');
disp('All required plots and details are generated.');
