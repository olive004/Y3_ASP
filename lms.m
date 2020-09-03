% ASP CW4 4.2.1 lms  
% filter order should be (Nw + 1), z & x should have shape ((Nw+1),1)
% output error: e               
% correction: mu e x
% input: x      
% desired: z      
% lr: mu      
% p: filt_order
function [y_hat, e_evol, w_evol] = lms(x, z, mu, filt_order)

    N = length(x);              % num of iterations
    perc_reduc_mu = 0.999;        % percent reduction in lr if the error is increasing
    threshold = 0.3;        % how much bigger the current error has to be compared to old
        
    w = zeros(filt_order, 1); 
    w_evol = zeros(filt_order, N);
    e_evol = zeros(N, 1);
    y_hat = zeros(N,1);

    for n= 1:N
        
        if n<=filt_order
            x_curr = x(1:filt_order);
        else
            x_curr = x(n-filt_order+1:n);
        end
        
        r_wx = w .* x_curr;
        y_hat(n) = sum(r_wx);
        
        % Theoretically: output_error(n) = desired(n) - y_hat(n);
        e = z(n) - y_hat(n);
    
        
        % GEAR SHIFT 4.3
        % Using rolling previous error instead of just e(n-1)
        window_length = 20;
        if n<=window_length
            past_e_window = window_length - (window_length - n);
        else
            past_e_window = window_length;
        end
        roll_past_e = sum(e_evol(n-past_e_window+1:n-1)) / past_e_window;
        is_significant =  (threshold * e) > roll_past_e;
        if is_significant
            if mu > 0.01 % Increasing this will cut the lr LESS, may mean fewer iterations
                mu = mu*threshold;
            else
                mu = mu*perc_reduc_mu;
            end
        end
        % Past input only 
%         if n>1
%             is_significant =  (threshold * e) > e_evol(n-1);
%             if is_significant
%                 if mu > 0.01
%                     mu = mu/10;
%                 else
%                     mu = mu*perc_reduc_mu;
%                 end
%             end
%         end
        
        
        
        
        
        % Correction: corr = lr * input * output_error(n);
        correction = mu * e * x_curr;
        w = w + correction; 
        
        w_evol(:, n) = w;
        e_evol(n) = e;
    
    end
    
    disp('Final lr mu: '); disp(mu);
end