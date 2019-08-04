clear all
close all
N = 35; % number of points
Fs = 1000; % sampling frequency
Fsig = 90; % signal frequency
columns = 7; %7
rows = 5;%5

t = (0:N-1)/N;
sig =sin(2*pi*Fsig*(0:N-1)/Fs);
subplot(5,1,1);
stem(t, sig); 
title('Input Signal');
ylabel('Magnitude');
%sig = 0:34;

tic
for v = 1:10
    X = DTFTfunction(sig);
end
tDTFT = toc;
display('Time Taken for DTFT')
disp(tDTFT)
subplot(5,1,2)
stem(20*log10(abs(X)), 'b');
ylabel('Magnitude (dB)');
title('DTFT');

tic
for v = 1:10
    % Reshape
    X = reshape(sig, [columns rows]).';
    % DTFT Columns
    for x = 1:columns
        X(:,x) = DTFTfunction(X(:,x));
    end

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
end

tFFT = toc;
display('Time Taken for FFT')
disp(tFFT);
% Display
hold on

subplot(5,1,3)
stem(20*log10(abs(abs(X))), 'r');
ylabel('Magnitude (dB)');
title('Cooley-Tukey FFT');
subplot(5,1,4)
stem(20*log10(abs(abs(fft(sig)))), 'g'); 
ylabel('Magnitude (dB)');
title('Matlab FFT');
%legend('DTFT', 'Cooley-Tukey', 'MATLAB fft()')


%reshape data
x = 1;
y = 1;
for i = 1:35
    DTFTGT(y, x) = sig(i);
    x = mod(x, columns)+1;
    y = mod(y, rows)+1;
end

for x = 1:7
        DTFTGT(:,x) = DTFTfunction(DTFTGT(:,x));
end

for y = 1:5
        DTFTGT(y,:) = DTFTfunction(DTFTGT(y,:));
end

for x = 1:7
    for y = 1:5
        k = mod(5*(x-1)+7*(y-1),35)+1;
        dataout(k) = DTFTGT(y,x);
    end
end
dataout;
hold on
subplot(5,1,5)
stem(20*log10(abs(dataout)))
ylabel('Magnitude (dB)');
title('Good-Thomas');