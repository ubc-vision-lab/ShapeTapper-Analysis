clear

out_path = 'D:\ShapeTapper-Analysis\';
% out_path = '';

patient_names = {'DF'};%, 'MC'};
% patient_names = {'MC1', 'MC2'};

conds = {'bounding_circle'};%{'in_shape', 'bounding_circle'};
metrics = {'amd_var'};%{'dplus_dminus', 'amd_var'};

num_trials = 100000;

statistics = {'task_correlation',...
              't_score_task','p_val_task','t_score_diff','p_val_diff',...
              't_diff_zscores','p_val_zscores','typeI_zscores',...
              't_diff_rsdt','p_val_rsdt','typeI_rsdt'};

for p=1:length(patient_names)
    
    for c=1:length(conds)
        
        for m=1:length(metrics)
            
            dat_fname = ['patient_and_control_scores_' conds{c} '_' metrics{m} '.mat'];
%             dat_fname = ['patient_and_control_scores_' conds{c} '_Simultaneous_2AFC_' metrics{m} '.mat'];
            
            try
                load(dat_fname);
            catch
                continue
            end
            
            n = length(subj_scores.Properties.RowNames) - 2;
            
            % task_labels = subj_scores.Properties.VariableNames(end);
            task_labels = subj_scores.Properties.VariableNames(end);
            dist_labels = subj_scores.Properties.VariableNames(1:end-1);  
            
            for t = 1:length(task_labels)

                t_diff = struct();
                for pp=1:length(patient_names)
                    t_diff.(patient_names{pp}) = cell2table(cell(length(statistics),length(dist_labels)), 'VariableNames', dist_labels, 'RowNames', statistics);
%                    t_diff.MC = cell2table(cell(length(statistics),length(dist_labels)), 'VariableNames', dist_labels, 'RowNames', statistics);
                end

                patient.task_data = subj_scores{patient_names(p),task_labels(t)} / 100;
                sample.task_data = subj_scores{3:end,task_labels(t)} / 100;

                [patient.task_t, patient.task_pval] = ttest_crawford(patient.task_data, sample.task_data);
                sample.task_ts = ttest_crawford(sample.task_data, sample.task_data);

                task_beta_pdf = fitdist(sample.task_data,'Beta');

                for d = 1:length(dist_labels)

                    fprintf(strjoin([patient_names{p},task_labels{t},dist_labels{d},"\n"],  ' '));

                    patient.dist_data = abs(subj_scores{patient_names(p),dist_labels(d)});
                    sample.dist_data = abs(subj_scores{3:end,dist_labels(d)});

                    [patient.dist_t, patient.dist_pval] = ttest_crawford(patient.dist_data, sample.dist_data);
                    sample.dist_ts = ttest_crawford(sample.dist_data, sample.dist_data);

                    if all(sample.dist_data < 1)
                        dist_beta_pdf = fitdist(sample.dist_data,'Beta');
                    else
                        dist_beta_pdf = fitdist(sample.dist_data,'Normal');
                    end

                    % Get t-score difference according to Crawford & Howell 1998
                    [patient.t_diff_zscores, patient.p_val_zscores] = tdiff_crawford(patient.task_data, sample.task_data, patient.dist_data, sample.dist_data);

                    % Get t-score according to the Revised Standardized Difference
                    % Test (RSDT), Crawford 2005a
                    [patient.t_diff_rsdt, patient.p_val_rsdt, patient.corr] = rsdt_crawford(patient.task_data, sample.task_data, patient.dist_data, sample.dist_data);

                    % Estimate Type I Error rates for t-diff for given beta PDFs
                    p_val_zscores = typeI_est_zscores(task_beta_pdf, dist_beta_pdf, n, num_trials);
                    p_val_rsdt = typeI_est_rsdt(task_beta_pdf, dist_beta_pdf, n, num_trials);

                    typeI_zscores = length(p_val_zscores(p_val_zscores < 2e-2)) / num_trials;
                    typeI_rsdt = length(p_val_rsdt(p_val_rsdt < 2e-2)) / num_trials;

                    % Save metrics to output table
                    t_diff.(patient_names{p}){'task_correlation', dist_labels{d}} = {patient.corr};

                    t_diff.(patient_names{p}){'t_score_task', dist_labels{d}} = {patient.task_t};
                    t_diff.(patient_names{p}){'p_val_task', dist_labels{d}} = {patient.task_pval};

                    t_diff.(patient_names{p}){'t_score_diff', dist_labels{d}} = {patient.dist_t};
                    t_diff.(patient_names{p}){'p_val_diff', dist_labels{d}} = {patient.dist_pval};

                    t_diff.(patient_names{p}){'t_diff_zscores', dist_labels{d}} = {patient.t_diff_zscores};
                    t_diff.(patient_names{p}){'p_val_zscores', dist_labels{d}} = {patient.p_val_zscores};
                    t_diff.(patient_names{p}){'typeI_zscores', dist_labels{d}} = {typeI_zscores};

                    t_diff.(patient_names{p}){'t_diff_rsdt', dist_labels{d}} = {patient.t_diff_rsdt};
                    t_diff.(patient_names{p}){'p_val_rsdt', dist_labels{d}} = {patient.p_val_rsdt};
                    t_diff.(patient_names{p}){'typeI_rsdt', dist_labels{d}} = {typeI_rsdt};


                end % distance measures

                fprintf(strjoin(["Writing output for",patient_names{p},task_labels{t},"\n"],  ' '));
                writetable(t_diff.(patient_names{p}), [out_path patient_names{p} '_' task_labels{t} '_' conds{c} '_' metrics{m} '_tdiff_analysis.xlsx'], 'WriteRowNames', true);

            end % tasks
        end % metrics
    end % conds
end % patients