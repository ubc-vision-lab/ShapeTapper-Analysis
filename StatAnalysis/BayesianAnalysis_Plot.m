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
% shapes = {'blake_01','blake_04','blake_06','blake_07',...
%           'blake_08','blake_10','blake_11','blake_12'}; 
shapes = {'blake_01','blake_03','blake_04','blake_06','blake_07',...
          'blake_08','blake_09','blake_10','blake_11','blake_12'};
num_shapes = length(shapes);

%%%%%%%%%%%%%%%% PLOT BY PATIENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=1:2
    for p=1:length(patients)

        out_path = ['C:\ShapeTapper-Analysis\bayesian_analysis\'];
        if ~exist(out_path, 'dir')
            mkdir(out_path);
        end

    %     bayes_fig = figure; hold on
    %     t = title([patients{p} ' - All Shapes']);
    %     set(t,'Interpreter', 'none')
        shapes_bf = array2table(zeros(num_shapes,6), ...
            'VariableNames', {'MedAxisVsNull',...
                              'CentVsNull', ...
                              'MedAxisVsCent', ...
                              'CentAndMedAxisVsNull', ...
                              'CentAndMedAxisVsCent', ...
                              'CentAndMedAxisVsMedaxis'},...
            'RowNames', shapes);
    %     shapes_bf = array2table(zeros(num_shapes,9), ...
    %         'VariableNames', {'MedAxis','MedAxisNoCent', 'Edge',...
    %                           'Cent', 'CentNoMedAxis', ...
    %                           'CentVsMedAxis', 'MedAxisVsCent',...
    %                           'CentNoMedAxisVsMedAxis', 'CentNoMedAxisVsMedAxisNoCent'},...
    %         'RowNames', shapes);
        px = 3;
        mm = px * 0.2482;

        for i=1:num_shapes

            in_path  = 'C:\ShapeTapper-Analysis\bayes_all_parts\';
%             '; %['C:\ShapeTapper-Analysis\' patients{p} '\bayesian_analysis\'];

            bayes_file = [in_path shapes{i} '_Patient_' patients{p} '_' num2str(t) '_bayes_rev.mat'];

            try
                data = load(bayes_file);
            catch
                continue
            end

    %         shapes_bf{i, 1} = data.medaxis_p / data.null_p;
    %         shapes_bf{i, 2} = data.centroid_p / data.null_p;
    %         shapes_bf{i, 3} = data.medaxis_p ./ data.centroid_p;
    %         shapes_bf{i, 4} = data.ma_plus_cent_p / data.null_p;
    %         shapes_bf{i, 5} = data.ma_plus_cent_p ./ data.centroid_p;
    %         shapes_bf{i, 6} = data.ma_plus_cent_p ./ data.medaxis_p;

    %         shapes_bf{i, 1} = data.medaxis_p(px) / data.null_p;
    %         shapes_bf{i, 2} = data.medaxis_solo_p(px) / data.null_p;
    %         shapes_bf{i, 3} = data.edge_p(px) / data.null_p;
    %         shapes_bf{i, 4} = data.centroid_p(px) / data.null_p;
    %         shapes_bf{i, 5} = data.centroid_solo_p(px) / data.null_p;
    %         shapes_bf{i, 6} = data.centroid_p(px) / data.medaxis_p(px);
    %         shapes_bf{i, 7} = data.medaxis_p(px) / data.centroid_p(px);
    %         shapes_bf{i, 8} = data.centroid_solo_p(px) / data.medaxis_p(px);
    %         shapes_bf{i, 9} = data.centroid_solo_p(px) / data.medaxis_solo_p(px);

    %         mdax_prob = data.medaxis_p / data.null_p;
    %         mdax_solo_prob = data.medaxis_solo_p / data.null_p;
    %         edge_prob = data.edge_p / data.null_p;
    %         cent_prob = data.centroid_p / data.null_p;

            mdax_prob = data.medaxis_p / data.null_p;
            cent_prob = data.centroid_p / data.null_p;
            mdax_plus_cent_prob = data.ma_plus_cent_p / data.null_p;

            shape_plot_name = erase(shapes{i},'_');
            out_path = [patients{p} '_rev\'];
            if ~exist(out_path, 'dir')
                mkdir(out_path);
            end

            xaxis_mm = (1:80)*0.2482;
            p1 = figure; hold on
            plot(xaxis_mm, ones(1, 80),'k', 'LineWidth', 2)
            plot(xaxis_mm, mdax_prob(1:80), 'r')
            plot(xaxis_mm, cent_prob(1:80), 'b');
            plot(xaxis_mm, mdax_plus_cent_prob(1:80), 'g')
            legend('Null','MedAxis Vs Null',...
                              'Cent Vs Null', ...
                              'MedAxis+Cent Vs Null');
            xlabel('Sigma (std of blur around MedAxis - NOT CENTROID) in mm');
            ylabel('Bayes Factor');
            title([patients{p} ' ' shape_plot_name ' task1 RefObjs vs Null - Bayes'])
            data.ma_plus_cent_p(data.ma_plus_cent_p==0) = NaN;
            data.medaxis_p(data.medaxis_p==0)  = NaN;
            data.centroid_p(data.centroid_p==0) = NaN;
            mdax_plus_cent_vs_solo_cent_prob = data.ma_plus_cent_p ./ data.centroid_p;
            mdax_plus_cent_vs_solo_mdax_prob = data.ma_plus_cent_p ./ data.medaxis_p;
            saveas(p1,[out_path patients{p} '_' shape_plot_name '_task' num2str(t) '_rev_bayes_factors_vs_null.png']);
            close(p1)

            p2=figure; hold on
            hold on
            xaxis_mm2 = (17:80)*0.2482;
            plot(xaxis_mm2, ones(1, 64),'k', 'LineWidth', 2)
            plot(xaxis_mm2, mdax_plus_cent_vs_solo_cent_prob(17:80), 'r')
            plot(xaxis_mm2, mdax_plus_cent_vs_solo_mdax_prob(17:80), 'b');
            legend('Null','MedAxis+Cent Vs Cent',...
                          'MedAxis+Cent Vs MedAxis');
            xlabel('Sigma (std of blur around MedAxis - NOT CENTROID) in mm');
            ylabel('Bayes Factor');
            title([patients{p} ' ' shape_plot_name ' task' num2str(t) ' MedAxis+Cent vs Solo RefObjs - Bayes'])
            saveas(p2,[out_path patients{p} '_' shape_plot_name '_task' num2str(t) '_rev_bayes_factors_sum_vs_ro.png']);
            close(p2)

    %         xaxis_mm = (1:60)*0.2482;
    %         
    %         plot(xaxis_mm, ones(1, 60),'k', 'LineWidth', 2)
    %         plot(xaxis_mm, edge_prob(1:60), 'r')
    %         plot(xaxis_mm, mdax_solo_prob(1:60), 'b');
    %         plot(cent_prob(1:60), 'g')
        end

    %     outname = strcat(patients{p},'_bayes_factors_1_',num2str(round(mm)),'mm.xlsx');

    %     writetable(shapes_bf,outname);
    %     f=get(gca,'Children');
    %     
    %     legend(f([1:3]),{'Medial Axis - Centroid',' All Edge','Null Hypothesis'});
    % 
    %     xticks(0:15)
    %     ylim([0 8])
    %     grid on
    %     set(gca, 'gridlinestyle', '--')
    %     xlabel("StdDev of Gauss smoothing in mm")
    %     ylabel("Bayes Factor (1 = NULL)")
    % 
    %     saveas(bayes_fig, strjoin([out_path, patients{p}, "_bayesian_probs.png"], ""));
    end

end

% %%%%%%%%%%%%%%%% PLOT BY SHAPE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% for i=1:num_shapes
% 
%     out_path = ['D:\ShapeTapper-Analysis\bayesian_analysis\'];
%     if ~exist(out_path, 'dir')
%         mkdir(out_path);
%     end
%     
%     bayes_fig = figure; hold on
%     t = title([shapes{i} ' - All Patients']);
%     set(t,'Interpreter', 'none')
%     
%     for p=1:length(patients)
%          
%         in_path  = ['D:\ShapeTapper-Analysis\' patients{p} '\bayesian_analysis\'];
%         
%         bayes_file = [in_path shapes{i} '_Patient_' patients{p} '_bayes.mat'];
% 
%         try
%             data = load(bayes_file);
%         catch
%             continue
%         end
% 
%         mdax_prob = data.medaxis_p / data.null_p;
%         mdax_solo_prob = data.medaxis_solo_p / data.null_p;
%         edge_prob = data.edge_p / data.null_p;
%         cent_prob = data.centroid_p / data.null_p;
%         
%         xaxis_mm = (1:60)*0.2482;
%         
%         plot(xaxis_mm, ones(1, 60),'k', 'LineWidth', 2)
%         plot(xaxis_mm, edge_prob(1:60), 'r')
%         
%         if strcmp(patients{p}, 'MC')
%             plot(xaxis_mm, mdax_solo_prob(1:60), 'g', 'LineWidth', 2);
%         elseif strcmp(patients{p}, 'DF')
%             plot(xaxis_mm, mdax_solo_prob(1:60), 'm', 'LineWidth', 2);
%         else
%             plot(xaxis_mm, mdax_solo_prob(1:60), 'b');
%         end
%         
% %         plot(cent_prob(1:60), 'g')
%     end
%     f=get(gca,'Children');
% 
%     if i == 2 % blake_03 has no data for DF
%         legend(f([1 7:9]),{'MC (Medial Axis - Centroid)','Control (Medial Axis - Centroid)','All Edge','Null Hypothesis'});
%     else
%         legend(f([1 4 7:9]),{'DF (Medial Axis - Centroid)', 'MC (Medial Axis - Centroid)','Control (Medial Axis - Centroid)','All Edge','Null Hypothesis'});
%     end
% 
%     xticks(0:15)
%     ylim([0 8])
%     grid on
%     set(gca, 'gridlinestyle', '--')
%     xlabel("StdDev of Gauss smoothing in mm")
%     ylabel("Bayes Factor (1 = NULL)")
% 
%     saveas(bayes_fig, strjoin([out_path, patients{p}, "_bayesian_probs.png"], ""));
% end