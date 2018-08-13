clear

% Start parallel pool (if parpool has started, simply ignore)
pool = gcp; 

% Get participant score data for t-tests
load('patient_and_control_scores_inshape.mat')

task_labels = subj_scores.Properties.VariableNames(end-4:end);
dist_labels = subj_scores.Properties.VariableNames(4:end-5);

patient_names = {'DF', 'MC'};

n = length(subj_scores.Properties.RowNames) - 2;
num_trials = 100000;

% output = struct();
% output.DF = cell2table(cell(length(task_labels),length(dist_labels)), 'VariableNames', dist_labels, 'RowNames', task_labels);
% output.MC = cell2table(cell(length(task_labels),length(dist_labels)), 'VariableNames', dist_labels, 'RowNames', task_labels);
% 
% tdiffs = struct();
% tdiffs.DF = cell2table(cell(length(task_labels),length(dist_labels)), 'VariableNames', dist_labels, 'RowNames', task_labels);
% tdiffs.MC = cell2table(cell(length(task_labels),length(dist_labels)), 'VariableNames', dist_labels, 'RowNames', task_labels);

tdiff_pvals_lq = struct();
tdiff_pvals_lq.DF = cell2table(cell(length(task_labels),length(dist_labels)), 'VariableNames', dist_labels, 'RowNames', task_labels);
tdiff_pvals_lq.MC = cell2table(cell(length(task_labels),length(dist_labels)), 'VariableNames', dist_labels, 'RowNames', task_labels);

tdiff_pvals_uq = struct();
tdiff_pvals_uq.DF = cell2table(cell(length(task_labels),length(dist_labels)), 'VariableNames', dist_labels, 'RowNames', task_labels);
tdiff_pvals_uq.MC = cell2table(cell(length(task_labels),length(dist_labels)), 'VariableNames', dist_labels, 'RowNames', task_labels);

% typeI_errs_twopct = struct();
% typeI_errs_twopct.DF = cell2table(cell(length(task_labels),length(dist_labels)), 'VariableNames', dist_labels, 'RowNames', task_labels);
% typeI_errs_twopct.MC = cell2table(cell(length(task_labels),length(dist_labels)), 'VariableNames', dist_labels, 'RowNames', task_labels);

for p = 1:2
    
    for t = 1:length(task_labels)
        
        patient.task_data = subj_scores{patient_names(p),task_labels(t)};
        sample.task_data = subj_scores{3:end,task_labels(t)};
        
        s_xx = sum((sample.task_data - mean(sample.task_data)).^2);
        
        for d = 1:length(dist_labels)
            
            fprintf(strjoin([patient_names{p},task_labels{t},dist_labels{d},"\n"],  ' '));
            
        	patient.dist_data = abs(subj_scores{patient_names(p),dist_labels(d)});
            sample.dist_data = abs(subj_scores{3:end,dist_labels(d)});
            
            s_yy = sum((sample.dist_data - mean(sample.dist_data)).^2);
            s_xy = sum((sample.dist_data - mean(sample.dist_data)).*(sample.task_data - mean(sample.task_data)));
            Tau = [s_xx, s_xy; s_xy, s_yy];
            [~, DI] = iwishrnd(Tau,n);
            
            control_means = [mean(sample.task_data); mean(sample.dist_data)];
            diff_mean_task = (patient.task_data - control_means(1));
            diff_mean_dist = (patient.dist_data - control_means(2));
            
            p_vals = zeros(1, num_trials);
%             [z_scores, p_vals] = deal(zeros(1, num_trials));
%             p_vals_all = zeros(n, num_trials);
%             z_scores = zeros(n, num_trials);
            parfor ii=1:num_trials
%                 fprintf([num2str(ii) '\n']);
                W_sample = iwishrnd(Tau, n, DI);
                T = chol(W_sample);
                zs = normrnd(0,1, [2,1]);
                stdevs_sample = control_means + (T*zs)/sqrt(n);

                z_score_task = diff_mean_task / W_sample(1,1);
                z_score_dist = diff_mean_dist / W_sample(2,2);
%                 z_score_task = (sample.task_data - control_means(1)) / W_sample(1,1);
%                 z_score_dist = (sample.dist_data - control_means(2)) / W_sample(2,2);

                rho_sample = W_sample(2,2) / sqrt(W_sample(1,1)*W_sample(2,2));

                z_score_sample = (z_score_task - z_score_dist) / sqrt(2 - 2*rho_sample);
%                 z_scores(ii) = z_score_sample;
                p_vals(ii) = normcdf(z_score_sample);
%                 p_vals_all(:,ii) = 1 - p_vals;
            end
            
            p_vals = sort(1 - p_vals);
            p_lq = p_vals(num_trials / 40);
            p_uq = p_vals(39 * num_trials / 40);
            
            tdiff_pvals_lq.(patient_names{p}){task_labels{t}, dist_labels{d}} = {p_lq};
            tdiff_pvals_uq.(patient_names{p}){task_labels{t}, dist_labels{d}} = {p_uq};
            
%             p_vals_means = mean(p_vals_all, 1);
%             p_vals_means = sort(p_vals_means);
%             p_lq = p_vals_means(num_trials / 40);
%             p_uq = p_vals_means(39 * num_trials / 40);
                                     
        end % distance measures
    end % tasks
end % patients