function sig = modelfun(coeffs,x)
%Model function for our motion

    %Lets say our motion consists of 10 frequencies 0.2Hz-4Hz
    freqs = linspace(0.2,5,length(coeffs)/2); % the coefficients should be an even number, since we have both 

    %Our first 10 coefficients are Weight_i  and last 10 are Phase_i
    %function is sig = sum(Weight_i*sin(2*pi*f*t+Phase_i));
    sig = zeros(size(x))';
    for i=1:numel(x)
        for j=1:numel(freqs)
            sig(i) = sig(i)+coeffs(j)*sin(2*pi*freqs(j)*x(i)+coeffs(j+numel(freqs)));
        end
    end
end