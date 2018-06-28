clear

% Load a.mat file with patient results
load('a.mat')

% Seperate the four columns into task categories 
% (a = 3 shape discrimination, b = 2 shape discrimination)
task_a = a(:,1);
task_b = nanmean(a(:,2:end),2);

pool = gcp; % Start parallel pool (if parpool has started, simply ignore)

%%%%% BETA Z-score %%%%%

% Find the beta distribution for both tasks
% Scale by 1/100 to convert values from percentage points to [0,1] range
pdf_a = fitdist(task_a(3:end)/100,'Beta');
pdf_b = fitdist(task_b(3:end)/100,'Beta');

% Calculate variance and standard deviation for task a
alpha_a = pdf_a.a;
beta_a  = pdf_a.b;

sample_size = length(task_a(3:end));
num_sets = 100000;
alphas_a = zeros(1,num_sets);
betas_a  = zeros(1,num_sets);
alphas_b = zeros(1,num_sets);
betas_b  = zeros(1,num_sets);
% beta_means_a  = zeros(1,num_sets);
% beta_means_b  = zeros(1,num_sets);
% beta_ameans_a  = zeros(1,num_sets);
% beta_ameans_b  = zeros(1,num_sets);
% beta_vars_a   = zeros(1,num_sets);
% beta_vars_b   = zeros(1,num_sets);
% beta_samplevars_a   = zeros(1,num_sets);
% beta_samplevars_b   = zeros(1,num_sets);
ppm1 = parfor_progressbar(num_sets,'Calculating alpha/beta values'); %create the progress bar
parfor i=1:num_sets
    ppm1.iterate(1); 
    rand_betas_a = random(pdf_a, [sample_size,1]);
    rand_betas_b = random(pdf_b, [sample_size,1]);
%     rand_betas_a = random(pdf_a, [50,1]);
%     rand_betas_b = random(pdf_b, [50,1]);
    
    rand_pdf_a = fitdist(rand_betas_a,'Beta');
    ralpha_a = rand_pdf_a.a;
    rbeta_a  = rand_pdf_a.b;
    alphas_a(i) = ralpha_a;
    betas_a(i)  = rbeta_a;
%     beta_means_a(i) = ralpha_a / (ralpha_a + rbeta_a);
%     beta_vars_a(i)  = (ralpha_a * rbeta_a) / ((ralpha_a + rbeta_a)^2 * (ralpha_a + rbeta_a + 1));
%     beta_samplevars_a(i) = var(rand_betas_a);
%     
    rand_pdf_b = fitdist(rand_betas_b,'Beta');
    alpha_b = rand_pdf_b.a;
    beta_b  = rand_pdf_b.b;
    alphas_b(i) = alpha_b;
    betas_b(i)  = beta_b;
%     beta_means_b(i) = alpha_b / (alpha_b + beta_b);
%     beta_vars_b(i)  = (alpha_b * beta_b) / ((alpha_b + beta_b)^2 * (alpha_b + beta_b + 1));
%     beta_samplevars_b(i) = var(rand_betas_b);
%     
%     beta_ameans_a(i)  = mean(rand_betas_a);
%     beta_ameans_b(i)  = mean(rand_betas_b);
end
close(ppm1);

[skewness_a, skewness_b, kurtosis_a, kurtosis_b] = deal(zeros(1, num_sets));
[skewnessr_a, skewnessr_b, kurtosisr_a, kurtosisr_b] = deal(zeros(1, num_sets));
ppm2 = parfor_progressbar(num_sets,'Calculating skewness/kurtosis values'); %create the progress bar
parfor jj = 1:num_sets
    ppm2.iterate(1); 
    fprintf([num2str(jj) '\n']);
    a = alphas_a(jj);
    b = betas_a(jj);
    skewness_a(jj) = (2*(b-a)*sqrt(a+b+1)) / ((a+b+2)*sqrt(a*b));
    kurtosis_a(jj) = ( (6*((a+b+1)*(a-b)^2 - a*b*(a+b+2))) / (a*b*(a+b+2)*(a+b+3)) ); 
    pdf_ab = makedist('Beta', 'a', a, 'b', b);
    rands_beta_ab = random(pdf_ab, [100000, sample_size]);
    skewnessr_a(jj) = mean(skewness(rands_beta_ab, 0, 2));
    kurtosisr_a(jj) = mean(kurtosis(rands_beta_ab, 0, 2));
   
    a = alphas_b(jj);
    b = betas_b(jj);
    skewness_b(jj) = (2*(b-a)*sqrt(a+b+1)) / ((a+b+2)*sqrt(a*b));
    kurtosis_b(jj) = ( (6*((a+b+1)*(a-b)^2 - a*b*(a+b+2))) / (a*b*(a+b+2)*(a+b+3)) ); 
    rands_beta_ab = random(pdf_ab, [100000, sample_size]);
    pdf_ab = makedist('Beta', 'a', a, 'b', b);
    skewnessr_b(jj) = mean(skewness(rands_beta_ab, 0, 2));
    kurtosisr_b(jj) = mean(kurtosis(rands_beta_ab, 0, 2));
end
close(ppm2);

skewness_explicit_task_a = skewness_a;
skewness_explicit_task_b = skewness_b;
skewness_montecarlo_task_a = skewnessr_a;
skewness_montecarlo_task_b = skewnessr_b;

kurtosis_explicit_task_a = kurtosis_a;
kurtosis_explicit_task_b = kurtosis_b;
kurtosis_montecarlo_task_a = kurtosisr_a;
kurtosis_montecarlo_task_b = kurtosisr_b;

save('skewness_kurtosis_1d_summary.mat',...
     'skewness_explicit_task_a', 'skewness_explicit_task_b',...
     'kurtosis_explicit_task_a', 'kurtosis_explicit_task_b',...
     'skewness_montecarlo_task_a', 'skewness_montecarlo_task_b',...
     'kurtosis_montecarlo_task_a', 'kurtosis_montecarlo_task_b',...
     'alphas_a', 'betas_a', 'alphas_b', 'betas_b')

n_steps_alpha = int32(100);
n_steps_beta  = int32(100);
% alpha_linspace = linspace(min(min(alphas_a),min(alphas_b)),max(max(alphas_a),max(alphas_b)), n_steps_alpha);
% beta_linspace  = linspace(min(min(betas_a),min(betas_b)),max(max(betas_a),max(betas_b)), n_steps_alpha);

alpha_min = min(min(alphas_a),min(alphas_b));
alpha_max = max(max(alphas_a),max(alphas_b));
alpha_diff = alpha_max - alpha_min;

beta_min = min(min(betas_a),min(betas_b));
beta_max = max(max(betas_a),max(betas_b));
beta_diff = beta_max - beta_min;

alpha_linspace = linspace(alpha_min, alpha_max, n_steps_alpha);
beta_linspace  = linspace(beta_min, beta_max, n_steps_beta);

skewness_2d = zeros(1, n_steps_alpha * n_steps_beta);
kurtosis_2d = zeros(1, n_steps_alpha * n_steps_beta);

ppm2 = parfor_progressbar(n_steps_alpha*n_steps_beta, 'Calculating skewness/kurtosis values');
parfor ii = 1:(n_steps_alpha*n_steps_beta)
    fprintf([num2str(ii) '\n']);ppm2.iterate(1);
%     if mod(ii,100)==0,  ppm2.iterate(100);  end
    alpha_idx = idivide(int32(ii)-1, n_steps_alpha) + 1;
    beta_idx  = mod(ii-1, n_steps_beta) + 1;
    alpha = double(alpha_idx)*alpha_diff/double(n_steps_alpha) + alpha_min;
    beta  = double(beta_idx)*beta_diff/double(n_steps_beta) + beta_min;
    pdf_ab = makedist('Beta', 'a', alpha, 'b', beta);
    rands_beta_ab = random(pdf_ab, [sample_size, num_sets]);
%     mean_rands_ab = mean(rands_beta_ab, 1);
    skewness_2d(ii) = mean(skewness(rands_beta_ab, 0, 1));
    kurtosis_2d(ii) = mean(kurtosis(rands_beta_ab, 0, 1)); 
end
close(ppm2);
% Reshape to 2d (alphas - columns, betas - rows)
skewness_2d_means = reshape(skewness_2d, [n_steps_alpha  n_steps_beta]);
kurtosis_2d_means = reshape(kurtosis_2d, [n_steps_alpha  n_steps_beta]);
% save('skewness_2d_random.mat', 'skewness_2d_random', 'alpha_linspace', 'beta_linspace')
% save('kurtosis_2d_random.mat', 'kurtosis_2d_random', 'alpha_linspace', 'beta_linspace')
% 
skewness_2d_explicit = zeros(n_steps_alpha, n_steps_beta);
kurtosis_2d_explicit = zeros(n_steps_alpha, n_steps_beta);

for j = 1:n_steps_alpha,
    fprintf([num2str(j) '\n']);
    for k = 1:n_steps_beta,
        a = alpha_linspace(j);
        b = beta_linspace(k);
        skewness_2d_explicit(j,k) = (2*(b-a)*sqrt(a+b+1)) / ((a+b+2)*sqrt(a*b));
        kurtosis_2d_explicit(j,k) = ( (6*((a+b+1)*(a-b)^2 - a*b*(a+b+2))) / (a*b*(a+b+2)*(a+b+3)) ); 
    end
end
skewness_2d_explicit = skewness_2d_explicit';
kurtosis_2d_explicit = kurtosis_2d_explicit';
% Reshape to 2d (alphas - columns, betas - rows)
% save('skewness_2d_explicit.mat', 'skewness_2d_explicit', 'alpha_linspace', 'beta_linspace')
% save('kurtosis_2d_explicit.mat', 'kurtosis_2d_explicit', 'alpha_linspace', 'beta_linspace')


alphas_task_a = alphas_a;
alphas_task_b = alphas_b;
betas_task_a  = betas_a;
betas_task_b  = betas_b;
save('skewness_2d_means.mat', 'skewness_2d_means', 'alpha_linspace', 'beta_linspace')
save('kurtosis_2d_means.mat', 'kurtosis_2d_means', 'alpha_linspace', 'beta_linspace')

save('skewness_kurtosis_2d_summary.mat',...
     'skewness_2d_means', 'kurtosis_2d_means',...
     'skewness_2d_explicit', 'kurtosis_2d_explicit',...
     'alpha_linspace', 'beta_linspace',...
     'alphas_task_a', 'alphas_task_b', 'betas_task_a', 'betas_task_b')

 save('skewness_kurtosis_2d_summary.mat',...
     'skewness_2d_explicit', 'kurtosis_2d_explicit',...
     'alpha_linspace', 'beta_linspace',...
     'alphas_task_a', 'alphas_task_b', 'betas_task_a', 'betas_task_b')

 
% 
% save('skewness_kurtosis_2d_summary.mat',...
%      'skewness_2d_random', 'skewness_2d_explicit',...
%      'kurtosis_2d_random', 'kurtosis_2d_explicit',...
%      'alpha_linspace', 'beta_linspace',...
%      'alphas_task_a', 'alphas_task_b', 'betas_task_a', 'betas_task_b')

% figure; hold on
% xs = linspace(0,1,100);
% pdf_mina_minb = pdf('Beta',xs, min(alphas_a), min(betas_a));
% plot(xs, pdf_mina_minb)
% pdf_mina_maxb = pdf('Beta',xs, min(alphas_a), max(betas_a));
% plot(xs, pdf_mina_maxb)
% pdf_maxa_minb = pdf('Beta',xs, max(alphas_a), min(betas_a));
% plot(xs, pdf_maxa_minb)
% pdf_maxa_maxb = pdf('Beta',xs, max(alphas_a), max(betas_a));
% plot(xs, pdf_maxa_maxb)
% 
% figure; hold on
% for i=1:10000
%     pdf_ab = pdf('Beta',xs, alphas_a(i), betas_a(i));
%     plot(xs, pdf_ab)
% end
% 
% nu = sample_size-1;
% [m,v] = tstat(nu);
% 
% z_score = (task_a(1) - m) / sqrt(v);
% 
% zscore(task_a)
% 
% mean_a = alpha_a / (alpha_a + beta_a);
% var_a = (alpha_a * beta_a) / ((alpha_a + beta_a)^2 * (alpha_a + beta_a + 1));
% std_a = sqrt(var_a);
% 
% % Calculate variance and standard deviation for task b
% alpha_b = pdf_b.a;
% beta_b  = pdf_b.b;
% mean_b = alpha_b / (alpha_b + beta_b);
% var_b = (alpha_b * beta_b) / ((alpha_b + beta_b)^2 * (alpha_b + beta_b + 1));
% std_b = sqrt(var_b);
% 
% % Calculate the zscore for DF and MC for tasks a and b
% % (Note the distribution is not necessarily normal!)
% zscore_beta_DF_a = (task_a(1)/100 - mean_a) / std_a; 
% zscore_beta_DF_b = (task_b(1)/100 - mean_b) / std_b; 
% 
% zscore_beta_MC_a = (task_a(2)/100 - mean_a) / std_a; 
% zscore_beta_MC_b = (task_b(2)/100 - mean_b) / std_b; 


% %%%%% HALF NORMAL Z-score %%%%%
% 
% % Get half normal distributions for tasks a and b
% % (We take 100 minus patient score to center scores on zero)
% halfnormal_a = fitdist(100-task_a(3:end),'HalfNormal');
% halfnormal_b = fitdist(100-task_b(3:end),'HalfNormal');
% 
% halfnormal_mean_a = halfnormal_a.mu + (halfnormal_a.sigma * (2/pi));
% halfnormal_var_a  = (halfnormal_a.sigma ^ 2) * (1 - (2/pi));
% halfnormal_std_a  = sqrt(halfnormal_var_a);
% 
% halfnormal_mean_b = halfnormal_b.mu + (halfnormal_b.sigma * (2/pi));
% halfnormal_var_b  = (halfnormal_b.sigma ^ 2) * (1 - (2/pi));
% halfnormal_std_b  = sqrt(halfnormal_var_b);
% 
% zscore_halfnormal_DF_a = (100-task_a(1) - halfnormal_mean_a) / halfnormal_std_a;
% zscore_halfnormal_DF_b = (100-task_a(1) - halfnormal_mean_b) / halfnormal_std_b;
% 
% zscore_halfnormal_MC_a = (100-task_a(2) - halfnormal_mean_a) / halfnormal_std_a;
% zscore_halfnormal_MC_b = (100-task_a(2) - halfnormal_mean_b) / halfnormal_std_b;
% 
% 
% %%%%% NORMAL Z-score %%%%%
% 
% % Get half normal distributions for tasks a and b
% % (We take 100 minus patient score to center scores on zero)
% normal_a = fitdist(task_a(3:end),'Normal');
% normal_b = fitdist(task_b(3:end),'Normal');
% 
% zscore_normal_DF_a = (task_a(1) - normal_a.mu) / normal_a.sigma;
% zscore_normal_DF_b = (task_b(1) - normal_b.mu) / normal_b.sigma;
% 
% zscore_normal_MC_a = (task_a(2) - normal_a.mu) / normal_a.sigma;
% zscore_normal_MC_b = (task_b(2) - normal_b.mu) / normal_b.sigma;
% 
% % Find skewness measure for normal distributions of tasks a and b
% skewness_a = skewness(task_a(3:end), 0);
% skewness_b = skewness(task_b(3:end), 0);

n_steps_alpha = int32(1000);
n_steps_beta  = int32(1000);

alpha_hist = histogram(alphas_task_a, n_steps_alpha);
as = alpha_hist.Values(:);
beta_hist = histogram(betas_task_a, n_steps_beta);
bs = beta_hist.Values(:)';
ab_hist = as * bs;
ab_hist = nthroot(ab_hist, 4);

figure;  subplot(1,2,1); grid on
s1 = surf(beta_linspace, alpha_linspace, kurtosis_2d_means, ab_hist);
s1.EdgeColor = 'none';
title({'Kurtosis of 100,000 sample means (N=29) by beta distribution parameters',...
       '(3 = normal distribution, >5 = moderate, >7 = severe)'})
xlabel('Beta values')
ylabel('Alpha values')
zlabel('Kurtosis of sample means')

subplot(1,2,2); grid on
s2 = surf(beta_linspace, alpha_linspace, skewness_2d_means, ab_hist);
s2.EdgeColor = 'none';
title({'Skewness of 100,000 sample means (N=29) by beta distribution parameters',...
        '(magnitude >0.31 = moderate, >0.7 = severe, >0.93 = very severe, >0.99 = extreme)'})
xlabel('Beta values')
ylabel('Alpha values')
zlabel('Skewness of sample means')


n_steps_alpha = int32(100);
n_steps_beta  = int32(100);

alpha_hist = histogram(alphas_task_a, n_steps_alpha);
as = alpha_hist.Values(:);
beta_hist = histogram(betas_task_a, n_steps_beta);
bs = beta_hist.Values(:)';
ab_hist = as * bs;
ab_hist = nthroot(ab_hist, 4);

figure;  subplot(2,2,1); grid on
s1 = surf(beta_linspace, alpha_linspace, kurtosis_2d_explicit + 3, ab_hist);
s1.EdgeColor = 'none';
title({'Mean kurtosis of 100,000 samples by beta distribution parameters (Explicit formula)',...
       '(3 = normal distribution, >5 = moderate, >7 = severe)'})
xlabel('Beta values')
ylabel('Alpha values')
zlabel('Kurtosis')
zlim([0 12])

subplot(2,2,2); grid on
s2 = surf(beta_linspace, alpha_linspace, skewness_2d_explicit, ab_hist);
s2.EdgeColor = 'none';
title({'Mean skewness of 100,000 samples by beta distribution parameters (Explicit formula)',...
        '(magnitude >0.31 = moderate, >0.7 = severe, >0.93 = very severe, >0.99 = extreme)'})
xlabel('Beta values')
ylabel('Alpha values')
zlabel('Skewness')
zlim([-3 1])

subplot(2,2,3); grid on
s1 = surf(beta_linspace, alpha_linspace, kurtosis_2d_random, ab_hist);
s1.EdgeColor = 'none';
title({'Mean kurtosis of 100,000 samples by beta distribution parameters (Monte Carlo method, N=29)',...
       '(3 = normal distribution, >5 = moderate, >7 = severe)'})
xlabel('Beta values')
ylabel('Alpha values')
zlabel('Kurtosis')
zlim([0 12])

subplot(2,2,4); grid on
s2 = surf(beta_linspace, alpha_linspace, skewness_2d_random, ab_hist);
s2.EdgeColor = 'none';
title({'Mean skewness of 100,000 samples by beta distribution parameters (Monte Carlo method, N=29)',...
        '(magnitude >0.31 = moderate, >0.7 = severe, >0.93 = very severe, >0.99 = extreme)'})
xlabel('Beta values')
ylabel('Alpha values')
zlabel('Skewness')
zlim([-3 1])

