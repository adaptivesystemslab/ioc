function [f, y, medVal, meanVal, intVal] = freqAnalysis(x, sampleFreq)
    % taking the FFT...
    [f, y] = takingFFT(x, sampleFreq);
    pdf = abs(y) .^ 2;
    
    medVal = medInt(f, pdf);
    meanVal = meanInt(f, pdf);
    intVal = trapz(pdf);
end

function [f, Yout] = takingFFT(emgSignal, freq)
	Fs = freq;                    % Sampling frequency
    T = 1/Fs;                     % Sample time
    L = length(emgSignal);        % Length of signal
    t = (0:L-1)*T;                % Time vector
    y = emgSignal;

    %plot(Fs*t(1:50),y(1:50))
    %%title('Signal Corrupted with Zero-Mean Random Noise')
    %xlabel('time (milliseconds)')

    NFFT = 2^nextpow2(L); % Next power of 2 from length of y
    Y = fft(y,NFFT)/L;
    f = Fs/2*linspace(0,1,NFFT/2+1);
    f = f';
    Yout = 2*abs(Y(1:NFFT/2+1));
    
    % Plot single-sided amplitude spectrum.
    figure;
    plot(f, Yout, 'r'); 
    title('Single-Sided Amplitude Spectrum of y(t)')
    xlabel('Frequency (Hz)')
    ylabel('|Y(f)|')

end

function medFreq = medInt(f, y)
    integralVal = cumtrapz(f, y);
    halfway = integralVal(end) / 2;
    
    intMatch = abs(integralVal - halfway);
    [C, I] = min(intMatch);
    medFreq = f(I);
end

function meanFreq = meanInt(f, y)
    num = f .* y;
    den = y;
    
    meanFreq = trapz(num) / trapz(den);
end