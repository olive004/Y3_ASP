% CW2.5.1 
% Average every step_size points together with input x scaled by alpha
% Input data should be in [data, samples] format
function x_smooth = smooth_ma(input, step_size, alpha)
    if size(input,1) == 1       % if only one sample supplied, reformat to [data,sample]
        input = input';
    end
    num_data = size(input,1);
    num_samples = size(input,2);
    excess = mod(num_data, step_size);
    
    % initialize output with 1/step_size length
    output_length = ceil(num_data/step_size);
    x_smooth = zeros(output_length, num_samples);
    
    % pad input with 0's
    input(end:end+(step_size-excess), :) = zeros((step_size-excess)+1, num_samples);

    % moving average w stepsize
    for sample_i = 1:num_samples
        for outp_i = 1:output_length-1
            inp_i = 1 + (outp_i-1)*step_size;
            x_smooth(outp_i, sample_i) = alpha * mean(input(inp_i:inp_i+step_size, sample_i));
            
        end
    end
end




