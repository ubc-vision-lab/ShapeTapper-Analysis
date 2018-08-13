clear

dat_path = "D:\ShapeTapper-Analysis\";
% dat_path = '';

shapes = {'blake_01','blake_03','blake_04','blake_06','blake_07',...
            'blake_08','blake_09','blake_10','blake_11','blake_12'};

conds = {'bounding_circle','in_shape'};

tasks = {};
% tasks = {'Simultaneous_2AFC'};%,'Sequential_2AFC','Oddball','N_Back'};

% % Parse D+/D- means by subject
% for c=1:length(conds)
%     for t=1:max(1,length(tasks))
%         
%         if isempty(tasks)
%            fname = strcat("Stats_allparts_shapemeans_",conds{c},"_dplus_dminus.xlsx");
%            outfile = strcat("patient_and_control_scores_",conds{c},"_dplus_dminus.mat");
%         else
%            fname = strcat("Stats_allparts_shapemeans_",conds{c},"_",tasks{t},"_dplus_dminus.xlsx");
%            outfile = strcat("patient_and_control_scores_",conds{c},"_",tasks{t},"_dplus_dminus.mat");
%         end
%            
%         [num,txt,raw] = xlsread(strcat(dat_path, fname));
%         num(num==0) = NaN;
%         
%         load(outfile)
%         load('a.mat', 'a')
% 
%         subjects = txt(2:end,1);
% 
%         dist_names = txt(1,2:end);
%         task_names = [dist_names {'Mean_Judgement_Accuracy'}];
% %         task_names = [strrep(a(1,:), '-', '_') {'Mean_2AFC', 'OutOfTwoAvg'}];
% %         task_names = [dist_names tasks];
% 
%         dist_data = num;
%         sample_task_data = nanmean(cell2mat(a(2:end,:)),2)/100;
%         patient_task_data = nanmean(subj_scores{1:2,8:11}, 2)/100;
% %         sample_task_data = [cell2mat(a(2:end,:)) mean(cell2mat(a(2:end,2:3)),2) mean(cell2mat(a(2:end,2:end)),2)];
% %         patient_task_data = [subj_scores{1:2, 10:13} mean(subj_scores{1:2, 11:12}, 2) subj_scores{1:2, 14}];
%         task_data = [patient_task_data; sample_task_data];
%         all_data = [dist_data task_data];
% 
%         clear subj_scores
%         subj_scores = array2table(all_data, 'VariableNames', task_names, 'RowNames', subjects');
% 
%         save(outfile, 'subj_scores')
%     end % tasks
% end % conds


clearvars -except dat_path conds tasks shapes

% Parse AMD and Variance by shape
for c=1:length(conds)
    for t=1:length(tasks)
        
        if isempty(tasks)
           fname = strcat("Stats_allparts_allshapes_",conds{c},"_amd_var.xlsx");
           outfile = strcat("patient_and_control_scores_",conds{c},"_amd_var.mat");
        else
           fname = strcat("Stats_allparts_allshapes_",conds{c},"_",tasks{t},"_amd_var.xlsx");
           outfile = strcat("patient_and_control_scores_",conds{c},"_",tasks{t},"_amd_var.mat");
        end
         
        [num,txt,raw] = xlsread(strcat(dat_path, fname));
        num(num==0) = NaN;
        
        load(outfile)
        load('a.mat', 'a')

        subjects = txt(2:end,1);

        idxs = [0:10:50];
        dist_names = txt(1,idxs+2);
        % Remove trailing _1 for shape labels
        for i=1:length(dist_names)
            under_idx = find(dist_names{i} == '_');
            dist_names{i}(under_idx(end):end) = [];
        end

        shape_tasks = {};
        for s=1:length(shapes)
            shape_tasks = [shape_tasks cellfun(@(x) strjoin({shapes{s},x},'_'),dist_names,'UniformOutput',false)];
        end

%         task_names = [strrep(a(1,:), '-', '_') {'Mean_2AFC', 'OutOfTwoAvg'}];
        task_names = [shape_tasks {'Mean_Judgement_Accuracy'}];

        dist_data = [];
        for ii=1:length(shapes)
            dist_data = [dist_data num(:,idxs+ii)];
        end

%         sample_task_data = cell2mat(a(2:end,3))/100;
%         patient_task_data = [0.617;0.55];
        sample_task_data = nanmean(cell2mat(a(2:end,:)),2)/100;
        patient_task_data = nanmean(subj_scores{1:2,61:64}, 2)/100;
%         sample_task_data = [cell2mat(a(2:end,:)) mean(cell2mat(a(2:end,2:3)),2) mean(cell2mat(a(2:end,2:end)),2)];
%         patient_task_data = [subj_scores{1:2, 10:13} mean(subj_scores{1:2, 11:12}, 2) subj_scores{1:2, 14}];
        task_data = [patient_task_data; sample_task_data];
        all_data = [dist_data task_data];

        clear subj_scores
        subj_scores = array2table(all_data, 'VariableNames', task_names, 'RowNames', subjects');

        save(outfile, 'subj_scores')
    end % tasks
end % conds