clear all


for k = 1:N
    X(k) = 0;
    for n = 0:N-1
        X(k) = X(k) + x(n+1)*exp(-1i * 2 * pi * k * n / N);
    end
end
stem(t, abs(real(X))); pause;
stem(t, abs(imag(X))); pause;