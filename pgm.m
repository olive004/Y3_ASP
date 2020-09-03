% CW3.0 script for periodogram
% Format: [samples, data]
function periodogram = pgm(x) 
    N = length(x);
    
    x_ft = fft(x);
    
    periodogram = (1/N) .* abs(x_ft).^2;

end