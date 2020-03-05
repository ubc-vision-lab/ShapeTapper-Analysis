clear; close all;

% patients = {'S01','S02','S03','S04','S05','S06','S07',...
%             'S08','S09','S10','S11','S12','S13','S14',...
%             'S15','S16','S17','S20','S21','S22','S23','S24',...
%             'S25','S26','S27','S28','S29','S30','S31','S32','S33',...
%             'S34','S35','S36','S37','S38','S39','S40'};
        
patients = {'S01','S02','S03','S04','S05','S06',...
    'S07','S08','S09','S10','S11','S12',...
    'S13','S14','S15','S16','S17','S18',...
    'S19','S20','S21','S22','S23','S24',...
    'S25','S26','S27','S28','S29','S30',...
    'S31','S32','S33','S34','S35','S36',...
    'S37','S38','S39','S40','S41','S42'};
        
%             'DF','MC'};
% patients = {'MC'};%{'DF','MC'};
% patients = {'S34', 'S35', 'S36', 'S37', 'S38', 'S39', 'S40',...
%             'S41', 'S42', 'S43', 'S44', 'S45', 'S46', 'S47',...
%             'S48', 'S49', 'S50', 'S51', 'S52'};

num_patients = length(patients);
shapes = {'blake_01','blake_04','blake_06','blake_07',...
          'blake_08','blake_10','blake_11','blake_12'}; 
% shapes = {'blake_01','blake_03','blake_04','blake_06','blake_07',...
%           'blake_08','blake_09','blake_10','blake_11','blake_12'};
num_shapes = length(shapes);

task1_bayes_unif_diffscores = nan(100000,num_shapes,num_patients);
task2_bayes_unif_diffscores = nan(100000,num_shapes,num_patients);

task1_bayes_obs_diffscores = nan(num_shapes,num_patients);
task1_bayes_obs_null_scores = nan(num_shapes,num_patients);

task2_bayes_obs_diffscores = nan(num_shapes,num_patients);
task2_bayes_obs_null_scores = nan(num_shapes,num_patients);


%%%%%%%%%%%%%%%% PLOT BY PATIENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=1:2
    for p=1:length(patients)

        out_path = ['C:\ShapeTapper-Analysis\bayesian_analysis\'];
        if ~exist(out_path, 'dir')
            mkdir(out_path);
        end

        for i=1:num_shapes

            in_path  = ['C:\ShapeTapper-Analysis\' patients{p} '\bayesian_analysis\'];
%             '; %'C:\ShapeTapper-Analysis\bayes_all_parts\';

            bayes_file = [in_path shapes{i} '_Patient_' patients{p} '_' num2str(t) '_bayes_unif_diff.mat'];

            try
                data = load(bayes_file);
            catch
                continue
            end

            if t==1
                task1_bayes_unif_diffscores(:,i,p) = data.centroid_p_unif - data.medaxis_p_unif;
                task1_bayes_obs_diffscores(i,p) = data.centroid_p_obs - data.medaxis_p_obs;
            elseif t==2
                task2_bayes_unif_diffscores(:,i,p) = data.centroid_p_unif - data.medaxis_p_unif;
                task2_bayes_obs_diffscores(i,p) = data.centroid_p_obs - data.medaxis_p_obs;
            end
            
        end % shapes
    end % patients

end % tasks