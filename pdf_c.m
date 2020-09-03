% Q1.3 Estimation of probability distributions
% OLD FUNCTION
% function [centers, scaled_counts] = pdf_c(x,bincount)  % where x is signal
%     
%     [counts, centers] = hist(x, bincount);
%     scaled_counts = counts/size(x,2);
% %     scaled_centers = 
% 
% %     figure;
% %     bar(centers, scaled_counts);
%     
% end

% copy of pdf.m
% Q1.3 gets and plots pdf; smoothing factor is by stepsize
function DF = pdf(x, n_bins, smoothing_factor)
    N = size(x, 2);
    pool_step = ceil(smoothing_factor);
    
    figure; hold on
    h = histogram(x,n_bins,'Normalization','pdf');
    heights = h.Values;
    centers = h.BinEdges(1:size(h.BinEdges,2)-1) + h.BinWidth/2;
    
    % Pool values within each pool_step; jumping average
    j = 1;          % preallocate
    heights_reduc = (pool_step):pool_step:length(heights);
    centers_reduc = (pool_step):pool_step:length(centers);
    for i=(pool_step):pool_step:length(heights)
        heights_reduc(j) = mean(heights((i-pool_step+1):i));
        centers_reduc(j) = mean(centers((i-pool_step+1):i));
        j=j+1;
        % Accounting for cut off end bit
        if (j==length(heights_reduc) & i<length(heights));
            heights_reduc = [heights_reduc mean(heights(i:length(heights)))];
            centers_reduc = [centers_reduc mean(centers(i:length(centers)))];
        end
    end

 
    n = length(centers_reduc);
    w = centers_reduc(2) - centers_reduc(1);
    t = linspace(centers_reduc(1)-w/2,centers_reduc(end)+w/2,n+1);    
    
    dt = diff(t);
    Fvals = cumsum([0,heights_reduc.*dt]); 
    
    F = spline(t, [heights_reduc(1), Fvals, heights_reduc(end)]); 
    DF = fnder(F);

    fnplt(DF, 'r', 2);
    hold off
    
end














