function [five_pct, two_pct, ts] = percent_typeones_diff(pdf, n, num_trials)
%percent_typeones_diff Percent Type I errors, differences in paired t-tests
% Note: as this funciton is vectorised use for small n 
% (requires [64 x n x num_trials] bytes)

% Generate pairs of random samples from Beta PDF
rands_beta_a = random(pdf, [num_trials, n+1, 2]);

% Calculate Crawford t-scores for each sample pair, get t-score differences
means = squeeze(mean(rands_beta_a(:, 1:end-1,:), 2));
stdevs = squeeze(std(rands_beta_a(:, 1:end-1,:), 0, 2));
t_factor = sqrt((n+1) / n);
ts_1 = (rands_beta_a(:,:,1) - means(:,1)) ./ (stdevs(:,1) * t_factor);
ts_2 = (rands_beta_a(:,:,2) - means(:,2)) ./ (stdevs(:,2) * t_factor);
t_diffs = ts_1 - ts_2;

% Perform Crawford t-test for each trial
stdev_t = std(t_diffs(:,1:end-1), 0, 2);
mean_t = mean(t_diffs(:,1:end-1), 2);
ts = (t_diffs(:,end) - mean_t) ./ (stdev_t * t_factor);

% Find percentage of p-vals under alpha range, i.e. num of Type I errors
five_pct = sum(abs(ts) > tinv(0.975, n)) / num_trials;
two_pct = sum(abs(ts) > tinv(0.99, n)) / num_trials;

end

