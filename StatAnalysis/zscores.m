clear

% Load a.mat file with patient results
load('a.mat')

% Seperate the four columns into task categories 
% (a = 3 shape discrimination, b = 2 shape discrimination)
task_a = a(:,1);
task_b = nanmean(a(:,2:end),2);


%%%%% BETA Z-score %%%%%

% Find the beta distribution for both tasks
% Scale by 1/100 to convert values from percentage points to [0,1] range
pdf_a = fitdist(task_a(3:end)/100,'Beta');
pdf_b = fitdist(task_b(3:end)/100,'Beta');

% Calculate variance and standard deviation for task a
alpha_a = pdf_a.a;
beta_a  = pdf_a.b;
mean_a = alpha_a / (alpha_a + beta_a);
var_a = (alpha_a * beta_a) / ((alpha_a + beta_a)^2 * (alpha_a + beta_a + 1));
std_a = sqrt(var_a);

% Calculate variance and standard deviation for task b
alpha_b = pdf_b.a;
beta_b  = pdf_b.b;
mean_b = alpha_b / (alpha_b + beta_b);
var_b = (alpha_b * beta_b) / ((alpha_b + beta_b)^2 * (alpha_b + beta_b + 1));
std_b = sqrt(var_b);

% Calculate the zscore for DF and MC for tasks a and b
% (Note the distribution is not necessarily normal!)
zscore_beta_DF_a = (task_a(1) - mean_a) / std_a; 
zscore_beta_DF_b = (task_b(1) - mean_b) / std_b; 

zscore_beta_MC_a = (task_a(2) - mean_a) / std_a; 
zscore_beta_MC_b = (task_b(2) - mean_b) / std_b; 


%%%%% HALF NORMAL Z-score %%%%%

% Get half normal distributions for tasks a and b
% (We take 100 minus patient score to center scores on zero)
halfnormal_a = fitdist(100-task_a(3:end),'HalfNormal');
halfnormal_b = fitdist(100-task_b(3:end),'HalfNormal');

halfnormal_mean_a = halfnormal_a.mu + (halfnormal_a.sigma * (2/pi));
halfnormal_var_a  = (halfnormal_a.sigma ^ 2) * (1 - (2/pi));
halfnormal_std_a  = sqrt(halfnormal_var_a);

halfnormal_mean_b = halfnormal_b.mu + (halfnormal_b.sigma * (2/pi));
halfnormal_var_b  = (halfnormal_b.sigma ^ 2) * (1 - (2/pi));
halfnormal_std_b  = sqrt(halfnormal_var_b);

zscore_halfnormal_DF_a = (100-task_a(1) - halfnormal_mean_a) / halfnormal_std_a;
zscore_halfnormal_DF_b = (100-task_a(1) - halfnormal_mean_b) / halfnormal_std_b;

zscore_halfnormal_MC_a = (100-task_a(2) - halfnormal_mean_a) / halfnormal_std_a;
zscore_halfnormal_MC_b = (100-task_a(2) - halfnormal_mean_b) / halfnormal_std_b;


%%%%% NORMAL Z-score %%%%%

% Get half normal distributions for tasks a and b
% (We take 100 minus patient score to center scores on zero)
normal_a = fitdist(task_a(3:end),'Normal');
normal_b = fitdist(task_b(3:end),'Normal');

zscore_normal_DF_a = (task_a(1) - normal_a.mu) / normal_a.sigma;
zscore_normal_DF_b = (task_b(1) - normal_b.mu) / normal_b.sigma;

zscore_normal_MC_a = (task_a(2) - normal_a.mu) / normal_a.sigma;
zscore_normal_MC_b = (task_b(2) - normal_b.mu) / normal_b.sigma;

% Find skewness measure for normal distributions of tasks a and b
skewness_a = skewness(task_a(3:end));
skewness_b = skewness(task_b(3:end));
