clear all;
close all;

% % -------------------- Parameters ----------------

[s , fs] = audioread("ECG.wav"); %loading ECG into an array and loading sampling frequency 
Ts = 1 / fs; %sample period
L = length(s); % length of the signal
s = transpose(s);
% creating time line for window function
k = 0:1:L-1;
blackman = 0.42 - 0.5*cos((2*pi*k)/(L-1)); 
hamming = 0.54-0.46*cos(2*pi*k/(L-1));
hanning = 0.5-0.5*cos(2*pi*k/(L-1));

x = randn(1,1000000)*sqrt(512);%generate white noise signal (normalise for default 512 point fft in pspectrum)

w = blackman; %window function
iirfc = 0.100; %cut off frequency for IIR
firfc = 50; %cut off frequency for FIR
FFT_size = fs; %FFT bin size
iirorder = 2; %IIR filter order.

m = 164; %number of FIR taps (N = 2m+1)
N = 2*m+1; %total number of FIR filter taps
adaptiveN = 5; %number of filter adaptive filter coefficients
mu = 0.001; %adaptive filter step size
ft = 50; %tonal interfiernece

%Task 1
% % -------------------- Time domain of the ECG.wav ----------------
% 
t = 0:Ts:(length(s)-1)*Ts; %generate discrete time values (nTs)
% 
% figure; plot(t,s); set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
% ylabel('Amplitude'); xlabel('Time [s]'); title('ECG.wav');
% grid on;
% 
% % % -------------------- Fast Fourier Transform (FFT) ----------------
% 
[magnitude, f] = power_spectrum(s, fs, FFT_size);

figure; plot(f , magnitude); grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14); xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]'); title('Spectrum of origianl ECG signal (own design)');
% 
% %%--------------------  pspectrum  ----------------
% figure;
% pspectrum(s, fs)               % spectrum of output signal (frequency response)
% ylabel('Gain (dB)');
% legend('ECG.wav');
% title('Spectrum of original ECG signal (pspectrum)');


%Task 2
% % % --------------------  Low frequency (drift) removal  ----------------
% 
% [T,B] = butter(iirorder,iirfc/(fs/2), 'high');
% disp(T);
% disp(B);
% 
% s_IIR = filter(T, B, s); 
% % figure; plot(t, s_IIR); set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
% % ylabel('Amplitude'); xlabel('Time [s]'); title('After MATLAB IIR ECG.wav');
% % grid on;
% % 
% % figure;
% % pspectrum(s_IIR, fs)               % spectrum of output signal (frequency response)
% % ylabel('Gain (dB)');
% % legend('FIR filter');
% % title('Spectrum of ECG after IIR filter (MATLAB)');
% 
% %--
% 
% [iir, b, a] = iir_filter(iirfc, fs, s);
% figure; plot(t, iir); set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
% ylabel('Amplitude'); xlabel('Time [s]'); title('After own desing IIR filter ECG.wav');
% grid on;
% 
% % %frequency response
% % figure; freqz(b,a,256,2000); title('Frequency Response of the own designed IIR filter');
% % 
% % figure; pspectrum(iir,fs)               % spectrum of output signal (frequency response)
% % ylabel('Gain (dB)');
% % legend('IIR filter');
% % title('Spectrum of ECG after IIR filter (Own designed');
% 
% % %_filter response
% % xf_IIR = filter(b, a , x);
% % figure; pspectrum(xf_IIR, fs)               % spectrum of output signal (frequency response)
% % ylabel('Gain (dB)');
% % legend('IIR filter');
% % title('Frequency Response IIR Filter');
% %--
% % 
% % % % --------------------  Tonal interference filtering  ----------------
% 
% for n = 1:m
%     h(n) = 2*(firfc/fs)*sin(n*2*pi*(firfc/fs))/(n*2*pi*(firfc/fs)); %calculate truncated impulse response for LP filter (+ve n)
% end
% h = [fliplr(h) 2*(firfc/fs) h]; %construct filter (add n = 0 coefficient for LP and -ve half)
% 
% w = hamming(N)';
% 
% hw = h.*w; %apply window to filter coefficients
% 
% %_________filter response
% xf_FIR = conv(hw,x);  %calculate filter output   
% figure; pspectrum(xf_FIR,fs)               % spectrum of output signal (frequency response)
% ylabel('Gain (dB)');
% legend('FIR filter');
% title('FIR Filter Frequency Response');
% %______
% 
% s_FIR = conv(hw,s_IIR);
% figure; plot(s_FIR); set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
% ylabel('Amplitude'); xlabel('Time [s]'); title('ECG.wav after FIR filter (conv)');
% grid on;
% 
% fir = fir_filter(firfc, fs, m, s_IIR);
% figure; plot(fir); set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
% ylabel('Amplitude'); xlabel('Time [s]'); title('ECG.wav after FIR filter (Own designed)');
% grid on;
% 
% figure;
% pspectrum(s_FIR, fs)               % spectrum of output signal (frequency response)
% ylabel('Gain (dB)');
% legend('FIR filter');
% title('Spectrum of ECG after FIR filter');
% 
% 
% 
% % % --------------------  Adaptive interference cancellation  ----------------
% e = adaptive_cancelation(fir, fs, ft, adaptiveN, mu);
% figure;plot(e)
% title('Estimated signal (error)');
% xlabel('no of iterations');
% ylabel('Amplitude');
% [p1,f] = pspectrum(s,fs);
% [p2,f1] = pspectrum(e,fs);
% 
% figure; plot(f,20*log(p1)); 
% hold on;
% plot(f, 20*log(p2));
% grid;
% title('Frequency Response');
% xlabel('Frequency [Hz]');
% ylabel('Amplitude [dB]');
% legend('ECG', 'ECG after filter');
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
% % --------------------   Multi-rate optimisation  ----------------

%%% graphs

% [p1,f] = pspectrum(s,fs);
% [p2,f1] = pspectrum(iir,fs);
% 
% figure; plot(f,20*log(p1)); 
% hold on;
% plot(f, 20*log(p2));
% grid;
% title('Frequency Response');
% xlabel('Frequency [Hz]');
% ylabel('Amplitude [dB]');
% 
% legend('ECG', 'ECG after IIR filter');
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
% 
% [p1,f] = pspectrum(s,fs);
% [p2,f1] = pspectrum(fir,fs);
% 
% figure; plot(f,20*log(p1)); 
% hold on;
% 
% plot(f, 20*log(p2));
% grid;
% title('Frequency Response');
% xlabel('Frequency [Hz]');
% ylabel('Amplitude [dB]');
% 
% legend('ECG', 'ECG after FIR filter');
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);

% %%adaptive filter compare
% [T,B] = butter(iirorder,iirfc/(fs/2), 'high');
% s_IIR = filter(T, B, s); 
% fir = fir_filter(firfc, fs, m, s);
