% ASP CW4: Sign LMS

function [x_hat, e_evol, w_evol] = sign_lms(x,mu,filt_order, is_sign_error, is_sign_regressor, is_gear_shift)
    % Make x a column vector (rows=1)
    if size(x, 2) == 1
        x = x';
    end
    N = length(x);
    % Make input delayed x sample (assuming a(1) == 1)
    x_decorr = [x; zeros(filt_order-1, N)];
    for order = 2:filt_order
        x_decorr(order,:) = [zeros(1,order), x(1:end-order)];
    end
    % Adaptive correction:
    perc_reduc_mu = 0.999;      % percent reduction in lr if the error is increasing
    threshold = 0.3;            % how much bigger the current error has to be compared to old
    
    
    w = zeros(filt_order, 1);
    w_evol = zeros(filt_order, N);
    e_evol = zeros(N, 1);
    x_hat = zeros(N,1);
    
    for n= 1:N
        x_curr = zeros(filt_order,1);
        for order = 1:filt_order
            if n<=filt_order
                x_curr(order) = x_decorr(order,n);
            else
                x_curr(order) = x_decorr(order, n - (filt_order - order));
            end
        end
        r_wx = w .* x_curr; 
        x_hat(n) = sum(r_wx);
    
        e = x(n) - x_hat(n);
        
        
        
        % GEAR SHIFT FROM OG LMS 
        % Using rolling previous error instead of just e(n-1)
        if is_gear_shift
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
        end
        
        
        
        
        
        % Correction: corr = lr * input * output_error(n);
        if is_sign_error && is_sign_regressor
            correction = mu * sign(e) * sign(x_curr);
        elseif is_sign_error && ~is_sign_regressor
            sign_x = sign_asp(x_curr);
            correction = mu * sign(e) * (x_curr);
        elseif ~is_sign_error && is_sign_regressor
            
            correction = mu * (e) * sign(x_curr);
        else
            correction = mu * e * x_curr;
        end
        w = w + correction; 

        w_evol(:, n) = w;
        e_evol(n) = e;
    end
    
    if is_gear_shift
        disp('Final lr mu: '); disp(mu);
    end
end