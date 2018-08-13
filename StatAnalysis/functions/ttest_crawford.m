function [tscore, pval] = ttest_crawford(patient_result, sample_data)
%ttest_crawford Performs t-test specified in Crawford 2006

% Perform Crawford t-test for each trial
n = length(sample_data);
stdev = nanstd(sample_data, 0); % use N-1 sample flag
tscore = (patient_result - nanmean(sample_data)) ./ (stdev * sqrt((n+1) / n));

% Get corresponding p-val for t-score
pval = tcdf(tscore,n);
end