% Load the audio file
[audioData, originalFs] = audioread('audioDSP.wav');

% Resample the audio to meet the Nyquist criterion
fc = 1500; % Cutoff frequency in Hz
fs = 2 * fc; % Sampling frequency in Hz
audioDataResampled = resample(audioData, fs, originalFs);

% Parameters for STFT and noise reduction
frameLength = 256; % Frame length for STFT
overlap = 128; % Overlap between frames
window = hamming(frameLength, 'periodic'); % Hamming window

% Compute the STFT of the resampled audio signal
[S, F, T] = stft(audioDataResampled, fs, 'Window', window, 'OverlapLength', overlap, 'FFTLength', frameLength);

% Estimate the noise power spectrum
noiseFrames = 1:10; % Use the first 10 frames to estimate noise
noiseEstimate = mean(abs(S(:, noiseFrames)).^2, 2);

% Perform spectral subtraction
S_magnitude = abs(S);
S_phase = angle(S);
S_clean_magnitude = max(S_magnitude - sqrt(noiseEstimate), 0); % Ensure no negative magnitudes

% Reconstruct the clean signal from magnitude and phase
S_clean = S_clean_magnitude .* exp(1i * S_phase);

% Perform inverse STFT
cleanAudioResampled = istft(S_clean, fs, 'Window', window, 'OverlapLength', overlap, 'FFTLength', frameLength);

% Ensure the reconstructed audio is real
cleanAudioResampled = real(cleanAudioResampled);

% Normalize the cleaned audio for saving
cleanAudioResampled = cleanAudioResampled / max(abs(cleanAudioResampled)); 

% Save the cleaned audio
audiowrite('cleaned_audio_fc1500.wav', cleanAudioResampled, fs);

% Plot the original and cleaned signals
t_original = (0:length(audioDataResampled)-1) / fs;
figure;
subplot(2, 1, 1);
plot(t_original, audioDataResampled);
title('Original Resampled Audio Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2, 1, 2);
t_clean = (0:length(cleanAudioResampled)-1) / fs;
plot(t_clean, cleanAudioResampled);
title('Cleaned Audio Signal');
xlabel('Time (s)');
ylabel('Amplitude');
