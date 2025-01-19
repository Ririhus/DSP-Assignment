
% Load the audio file
[input_signal, fs] = audioread('audioDSP.wav');

% Normalize the input signal
input_signal = input_signal / max(abs(input_signal));

% Perform FFT to get the frequency spectrum
N = length(input_signal);
frequencies = (0:N-1) * (fs / N); % Frequency axis
spectrum = abs(fft(input_signal)); % Magnitude spectrum

% Plot the frequency spectrum
figure;
plot(frequencies(1:floor(N/2)), spectrum(1:floor(N/2))); 
% Single-sided spectrum
title('Frequency Spectrum of Input Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% Compute and plot the spectrogram
figure;
spectrogram(input_signal, 256, 128, 512, fs, 'yaxis');
title('Spectrogram of Input Signal');
colormap jet;
