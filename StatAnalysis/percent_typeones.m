function [five_pct, two_pct] = percent_typeones(pdf, n, num_trials)
%percent_typeones Percentage of Type I errors using Crawfords t-test

% Generate random samples from Beta PDF
rands = random(pdf, [num_trials, n+1]);

% Perform Crawford t-test for each trial
stdevs = std(rands, 0, 2);
means = mean(rands(:,1:end-1), 2);
ts = (rands(:,end) - means) ./ (stdevs * sqrt((n+1) / n));

% Find percentage of p-vals under alpha range, i.e. num of Type I errors
five_pct = sum(abs(ts) > tinv(0.975, n)) / num_trials;
two_pct = sum(abs(ts) > tinv(0.99, n)) / num_trials;

end

