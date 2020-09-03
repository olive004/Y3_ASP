% Function for getting ave PSD via divisions
function [psd_ave] = get_ave_psd(x, window_size)
    N = length(x);
    psd_wind = zeros(window_size, N/window_size);
    
    for div = 1:N/window_size
        i1 = (div-1)*window_size+1;
        i2 = div*window_size;
        
        psd_wind(1:window_size, div) = pgm(x(i1:i2));

    end
    psd_ave = mean(psd_wind,2);


end




