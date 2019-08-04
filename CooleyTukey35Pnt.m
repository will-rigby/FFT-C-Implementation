clear all
close all
N = 35; % number of points
Fs = 1000; % sampling frequency
Fsignal = 90; % signalnal frequency
columns = 5; %7
rows = 7;%5

t = (0:N-1)/N;
signal =sin(2*pi*Fsignal*(0:N-1)/Fs)
subplot(4,1,1);
stem(t, signal); 
title('Input signalnal');
ylabel('Magnitude');
%signal = 0:34;

count = 0;

tic
%for v = 1:100
    
    X = DTFTfunction(signal);
%end
tDTFT = toc;
display('Time Taken for DTFT')
disp(tDTFT)
subplot(4,1,2);
stem(20*log10(abs(X)), 'b');
ylabel('Magnitude (dB)');
title('DTFT');

tic
%for v = 1:100
    % Reshape
    X = reshape(signal, [columns rows]).'
    % DTFT Columns
    for x = 1:columns
        X(:,x) = DTFTfunction(X(:,x));
    end
    X

    % Twiddle Factor Operations
    for x = 1:columns
        for y = 1:rows
            X(y,x)=X(y,x)*exp((-2*pi*(y-1)*(x-1)/N)*1i);
        end
    end

    %DTFT Rows
    for y = 1:rows
        X(y,:) = DTFTfunction(X(y,:));
    end

    % Reshape
    X = reshape(X, 1, N);
%end

tFFT = toc;
display('Time Taken for FFT')
disp(tFFT);
% Display
hold on

subplot(4,1,3);
stem(20*log10(abs(abs(X))), 'r');
ylabel('Magnitude (dB)');
title('Cooley-Tukey FFT');
subplot(4,1,4);
stem(20*log10(abs(abs(fft(signal)))), 'g'); 
ylabel('Magnitude (dB)');
title('Matlab FFT');
%legend('DTFT', 'Cooley-Tukey', 'MATLAB fft()')

