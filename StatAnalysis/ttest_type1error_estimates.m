clear

pool = gcp; % Start parallel pool (if parpool has started, simply ignore)

% Load alpha/beta values fitting Task1_Accuracy
load('skewness_kurtosis_1d_summary.mat')
clearvars -except alphas_a alphas_b betas_a betas_b

% Load mat file with patient results
load('patient_and_control_scores_inshape.mat')

% Seperate the columns into task categories 
task_a = subj_scores.Task1_Accuracy;
% task_b = subj_scores.MedialAxis_AMD_Observed;

% Get sample size (-2 for MC and DF) and Monte Carlo simulation count
n = length(task_a) - 2;
num_trials = 1000000;

% Set up measurements for alpha-beta linspaces
n_steps_alpha = int32(50);
n_steps_beta  = int32(50);
alpha_min = min(min(alphas_a),min(alphas_b));
alpha_max = max(max(alphas_a),max(alphas_b));
alpha_diff = alpha_max - alpha_min;
beta_min = min(min(betas_a),min(betas_b));
beta_max = max(max(betas_a),max(betas_b));
beta_diff = beta_max - beta_min;

% Create alpha/beta probability coloring for surface plots
alpha_ax = linspace(alpha_min, alpha_max, n_steps_alpha);
beta_ax = linspace(beta_min, beta_max, n_steps_beta);

alpha_hist = histogram(alphas_a, n_steps_alpha);
as = alpha_hist.Values(:);
beta_hist = histogram(betas_a, n_steps_beta);
bs = beta_hist.Values(:)';
ab_hist = as * bs;
% ab_hist = nthroot(ab_hist, 4);
% ab_hist = ab_hist / sum(sum(ab_hist));


%% Calculate and plot Type I error estimates for Crawford t-scores

[fivepct, twopct] = deal(zeros(1, n_steps_alpha * n_steps_beta));
parfor ii=1:(n_steps_alpha * n_steps_beta)
    fprintf([num2str(ii) '\n']);
    
    alpha_idx = idivide(int32(ii)-1, n_steps_alpha) + 1;
    beta_idx  = mod(ii-1, n_steps_beta) + 1;
    alpha = double(alpha_idx)*alpha_diff/double(n_steps_alpha) + alpha_min;
    beta  = double(beta_idx)*beta_diff/double(n_steps_beta) + beta_min;
    
    pdf_beta = makedist('Beta', 'a', alpha, 'b', beta);
    [fivepct(ii), twopct(ii)] = percent_typeones(pdf_beta, n, num_trials);
end

typeonepct_five = reshape(fivepct, [n_steps_alpha  n_steps_beta]);
typeonepct_two = reshape(twopct, [n_steps_alpha  n_steps_beta]);

figure;  hold on; grid on
s5 = surf(beta_ax, alpha_ax, typeonepct_five * 100, ab_hist);
s5.EdgeColor = 'none';
s2 = surf(beta_ax, alpha_ax, typeonepct_two * 100, ab_hist);
s2.EdgeColor = 'none';
title({'Type I error rate in Crawford t-test for beta distributions';...
       '1,000,000 Monte Carlo simulations (N=29) for each alpha/beta pair'})
xlabel('Beta values')
ylabel('Alpha values')
zlabel('Percentage of Type I Errors')
ztickformat('percentage')
legend('p < 0.05 (above)', 'p < 0.02 (below)')


%% Calculate and plot Type I error estimates for Crawford difference t-scores

[fivepct, twopct] = deal(zeros(1, n_steps_alpha * n_steps_beta));
parfor ii=1:(n_steps_alpha * n_steps_beta)
    fprintf([num2str(ii) '\n']);
    
    alpha_idx = idivide(int32(ii)-1, n_steps_alpha) + 1;
    beta_idx  = mod(ii-1, n_steps_beta) + 1;
    alpha = double(alpha_idx)*alpha_diff/double(n_steps_alpha) + alpha_min;
    beta  = double(beta_idx)*beta_diff/double(n_steps_beta) + beta_min;
    
    pdf_beta = makedist('Beta', 'a', alpha, 'b', beta);
    [fivepct(ii), twopct(ii)] = percent_typeones_diff(pdf_beta, n, num_trials);
end

typeonepct_five = reshape(fivepct, [n_steps_alpha  n_steps_beta]);
typeonepct_two = reshape(twopct, [n_steps_alpha  n_steps_beta]);

figure;  hold on; grid on
s5 = surf(beta_ax, alpha_ax, typeonepct_five * 100, ab_hist);
s5.EdgeColor = 'none';
s2 = surf(beta_ax, alpha_ax, typeonepct_two * 100, ab_hist);
s2.EdgeColor = 'none';
title({'Type I error rate in Crawford t-test difference scores for beta dists';...
       '1,000,000 Monte Carlo simulations (N=29) for each alpha/beta pair'})
xlabel('Beta values')
ylabel('Alpha values')
zlabel('Percentage of Type I Errors')
ztickformat('percentage')
legend('p < 0.05 (above)', 'p < 0.02 (below)')


%% Calculate and plot Type I error estimates for Crawford difference t-scores by N

alpha_mean = mean(alphas_a);
beta_mean  = mean(betas_a);

n_range = 100;
[fivepct_ns, twopct_ns] = deal(zeros(1, n_range));
pdf_beta_mean = makedist('Beta', 'a', alpha_mean, 'b', beta_mean);

for nn=1:n_range
    fprintf([num2str(nn) '\n']);
    [fivepct_ns(nn), twopct_ns(nn)] = percent_typeones_diff(pdf_beta_mean, nn, num_trials);
end

fivepct_ns(1) = 100;
twopct_ns(1) = 100;

figure;  hold on; grid on
plot(1:n_range, fivepct_ns * 100, 'LineWidth', 3);
plot(1:n_range, twopct_ns * 100, 'LineWidth', 3);
title({'Type I error rate in Crawford t-test difference scores for beta dists (mean \alpha,\beta)';...
       '1,000,000 Monte Carlo simulations for varying N [1 .. 100]'})
xlim([1 80])
xticks(0:5:80)
ylim([0 10])
xlabel('N')
ylabel('Percentage of Type I Errors')
ytickformat('percentage')
legend('p < 0.05', 'p < 0.02')


%% Calculate and plot histogram differences Crawford t-scores

[~,~,t_diffs] = percent_typeones_diff(pdf_beta_mean, n, num_trials);
fivepct_t = prctile(t_diffs,97.5);
twopct_t  = prctile(t_diffs,99);

figure;  hold on; grid on
hist(t_diffs, 100);
ymax = ylim;
plot([twopct_t twopct_t], [0 ymax(2)], 'r', 'LineWidth', 2);
plot([-twopct_t -twopct_t], [0 ymax(2)], 'r', 'LineWidth', 2);
plot([fivepct_t fivepct_t], [0 ymax(2)], 'g', 'LineWidth', 2);
plot([-fivepct_t -fivepct_t], [0 ymax(2)], 'g', 'LineWidth', 2);
title({'Distribution of Crawford t-test difference scores for beta dists (mean \alpha,\beta)';...
       '1,000,000 Monte Carlo simulations (N=26)'})
xlabel('t-difference')
xticks(-14:2:14)
ylabel('Frequency / 1000000 trials')
f=get(gca,'Children');
legend(f([1 3]),{'p < 0.05','p < 0.02'});