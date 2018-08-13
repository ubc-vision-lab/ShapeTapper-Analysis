clear

patients = {'DF','MC',...
            'S01','S02','S03','S04','S05','S06',...
            'S07','S08','S09','S10','S11','S12',...
            'S13','S14','S15','S16','S17','S18',...
            'S19','S20','S21','S22','S23','S24',...
            'S25','S26','S27','S28','S29','S30',...
            'S31','S32','S33','S34','S35','S36',... 
            'S37','S38','S39','S40','S41','S42',...
            'S43','S44','S45','S46','S47','S48',...
            'S49','S50','S51','S52','S53','S54'};

% patients = {'DF'};
num_patients = length(patients);

analysis_conds = {'in_shape','bounding_circle'};

tasks = {};
% tasks = {'Simultaneous_2AFC','Sequential_2AFC','Oddball','N_Back'};

shapes = {'blake_01','blake_03','blake_04','blake_06','blake_07',...
          'blake_08','blake_09','blake_10','blake_11','blake_12'};
      
num_shapes = length(shapes);

% Allocate structs for all subjs summary
all_ma = struct();
all_edge = struct();
all_cent = struct();

% Init distance structs for all subjs
[all_ma.dist_std, all_ma.dist_avg] = deal(zeros(num_patients,num_shapes));
[all_edge.dist_std, all_edge.dist_avg] = deal(zeros(num_patients,num_shapes));
[all_cent.dist_std, all_cent.dist_avg] = deal(zeros(num_patients,num_shapes));

[all_ma.obs_dist_var, all_ma.obs_dist_avg] = deal(zeros(num_patients,num_shapes));
[all_edge.obs_dist_var, all_edge.obs_dist_avg] = deal(zeros(num_patients,num_shapes));
[all_cent.obs_dist_var, all_cent.obs_dist_avg] = deal(zeros(num_patients,num_shapes));

% Init Dplus/Dminus structs for all subjs
[all_ma.dplus_std, all_ma.dplus_avg, all_ma.dminus_std, all_ma.dminus_avg] = deal(zeros(num_patients,num_shapes));
[all_edge.dplus_std, all_edge.dplus_avg, all_edge.dminus_std, all_edge.dminus_avg] = deal(zeros(num_patients,num_shapes));
[all_cent.dplus_std, all_cent.dplus_avg, all_cent.dminus_std, all_cent.dminus_avg] = deal(zeros(num_patients,num_shapes));

[all_ma.dmaxdev_avg, all_edge.dmaxdev_avg, all_cent.dmaxdev_avg] = deal(zeros(num_patients,num_shapes));
[all_ma.dmaxdev_std, all_edge.dmaxdev_std, all_cent.dmaxdev_std] = deal(zeros(num_patients,num_shapes));

[all_ma.obs_dplus, all_ma.obs_dminus] = deal(zeros(num_patients,num_shapes));
[all_edge.obs_dplus, all_edge.obs_dminus] = deal(zeros(num_patients,num_shapes));
[all_cent.obs_dplus, all_cent.obs_dminus] = deal(zeros(num_patients,num_shapes));


for c=1:length(analysis_conds)
    
    for t=1:max(length(tasks),1)
        
        for p=1:length(patients)
        
            if isempty(tasks)
                fprintf('Starting %s, %s, %s', patients{p}, analysis_conds{c});
                in_path  = ['D:\ShapeTapper-Analysis\' patients{p} '\spatial_analysis\' analysis_conds{c} '\'];
            else
                fprintf('Starting %s, %s', patients{p}, analysis_conds{c}, tasks{t});
                in_path  = ['D:\ShapeTapper-Analysis\' patients{p} '\spatial_analysis\' analysis_conds{c} '\' tasks{t} '\'];
            end
            
            out_path = ['D:\ShapeTapper-Analysis\' patients{p} '\spatial_analysis\'];

            if ~exist(out_path, 'dir')
                mkdir(out_path);
            end

            n = zeros(num_shapes,1);

            [ma_dist_std, ma_dist_avg] = deal(zeros(num_shapes,1));
            [edge_dist_std, edge_dist_avg] = deal(zeros(num_shapes,1));
            [cent_dist_std, cent_dist_avg] = deal(zeros(num_shapes,1));

            [ma_dist_obs, edge_dist_obs, cent_dist_obs] = deal(zeros(num_shapes,1));
            [ma_dist_obs_var, edge_dist_obs_var, cent_dist_obs_var] = deal(zeros(num_shapes,1));

            [medaxis_var, medaxis_amd, edge_var, edge_amd, centroid_var, centroid_amd] = deal(zeros(num_shapes,1));
            [medaxis_dplus, medaxis_dminus, edge_dplus, edge_dminus, centroid_dplus, centroid_dminus] = deal(zeros(num_shapes,1));
            [medaxis_result, edge_result, centroid_result] = deal(cell(num_shapes,1));

            [ma_dplus_std, ma_dplus_avg, ma_dminus_std, ma_dminus_avg] = deal(zeros(num_shapes,1));
            [edge_dplus_std, edge_dplus_avg, edge_dminus_std, edge_dminus_avg] = deal(zeros(num_shapes,1));
            [cent_dplus_std, cent_dplus_avg, cent_dminus_std, cent_dminus_avg] = deal(zeros(num_shapes,1));

            [ma_dmaxdev_avg, edge_dmaxdev_avg, cent_dmaxdev_avg] = deal(zeros(num_shapes,1));
            [ma_dmaxdev_std, edge_dmaxdev_std, cent_dmaxdev_std] = deal(zeros(num_shapes,1));

            [ma_dplus_obs, edge_dplus_obs, cent_dplus_obs] = deal(zeros(num_shapes,1));
            [ma_dminus_obs, edge_dminus_obs, cent_dminus_obs] = deal(zeros(num_shapes,1));

            for i=1:num_shapes
                fprintf('.');
                
                if isempty(tasks)
                    spat_file = [in_path shapes{i} '_Patient_' patients{p} '_spatial_analysis_' analysis_conds{c} '.mat'];
                else 
                    spat_file = [in_path shapes{i} '_Patient_' patients{p} '_spatial_analysis_' analysis_conds{c} '_' tasks{t} '.mat'];
                end
                
                try
                    data = load(spat_file);
                catch
                    continue
                end

                n(i) = data.n_points;

%                 medaxis_var(i)  = KernelDistApproximator(data.uniform_medaxis_data(1,:)', data.observed_medaxis_data(1)); % variance
%                 medaxis_amd(i)  = KernelDistApproximator(data.uniform_medaxis_data(3,:)', data.observed_medaxis_data(3)); % amd 
%                 edge_var(i)     = KernelDistApproximator(data.uniform_edge_data(1,:)', data.observed_edge_data(1)); % variance
%                 edge_amd(i)     = KernelDistApproximator(data.uniform_edge_data(3,:)', data.observed_edge_data(3)); % amd    
%                 centroid_var(i) = KernelDistApproximator(data.uniform_centroid_data(1,:)', data.observed_centroid_data(1)); % variance
%                 centroid_amd(i) = KernelDistApproximator(data.uniform_centroid_data(3,:)', data.observed_centroid_data(3)); % amd    

                ma_dist_std(i)   = std(data.uniform_medaxis_data(3,:),[],2);
                ma_dist_avg(i)   = mean(data.uniform_medaxis_data(3,:),2);
                edge_dist_std(i) = std(data.uniform_edge_data(3,:),[],2); 
                edge_dist_avg(i) = mean(data.uniform_edge_data(3,:),2);
                cent_dist_std(i) = std(data.uniform_centroid_data(3,:),[],2); 
                cent_dist_avg(i) = mean(data.uniform_centroid_data(3,:),2);

                if isempty(data.observed_medaxis_data)
                    ma_dist_obs(i) = NaN;
                    ma_dist_obs_var(i) = NaN;
                else
                    ma_dist_obs(i) = data.observed_medaxis_data(3);
                    ma_dist_obs_var(i) = data.observed_medaxis_data(1);
                end
                
                if isempty(data.observed_edge_data)
                    edge_dist_obs(i) = NaN;
                    edge_dist_obs_var(i) = NaN;
                else
                    edge_dist_obs(i) = data.observed_edge_data(3);
                    edge_dist_obs_var(i) = data.observed_edge_data(1);
                end
                
                if isempty(data.observed_centroid_data)
                    cent_dist_obs(i) = NaN;
                    cent_dist_obs_var(i) = NaN;
                else
                    cent_dist_obs(i) = data.observed_centroid_data(3);
                    cent_dist_obs_var(i) = data.observed_centroid_data(1);
                end

%                 [obs_ma, gen_ma]     = parse_cdf_struct(data.medaxis_cdf);
%                 [obs_edge, gen_edge] = parse_cdf_struct(data.edge_cdf);
%                 [obs_cent, gen_cent] = parse_cdf_struct(data.centroid_cdf);
% 
%                 ma_dplus_obs(i)    = obs_ma.dplus;
%                 ma_dminus_obs(i)   = obs_ma.dminus;
%                 edge_dplus_obs(i)  = obs_edge.dplus;
%                 edge_dminus_obs(i) = obs_edge.dminus;
%                 cent_dplus_obs(i)  = obs_cent.dplus;
%                 cent_dminus_obs(i) = obs_cent.dminus;

                if data.medaxis_cdf_obs == [0,0]
                    data.medaxis_cdf_obs = [NaN,NaN];
                end
                if data.edge_cdf_obs == [0,0]
                    data.edge_cdf_obs = [NaN,NaN];
                end
                if data.centroid_cdf_obs == [0,0]
                    data.centroid_cdf_obs = [NaN,NaN];
                end
                ma_dplus_obs(i)    = data.medaxis_cdf_obs(1);
                ma_dminus_obs(i)   = data.medaxis_cdf_obs(2);
                edge_dplus_obs(i)  = data.edge_cdf_obs(1);
                edge_dminus_obs(i) = data.edge_cdf_obs(2);
                cent_dplus_obs(i)  = data.centroid_cdf_obs(1);
                cent_dminus_obs(i) = data.centroid_cdf_obs(2);

%                 ma_dmaxdev_avg(i)   = gen_ma.dmaxdev_avg;
%                 edge_dmaxdev_avg(i) = gen_edge.dmaxdev_avg;
%                 cent_dmaxdev_avg(i) = gen_cent.dmaxdev_avg;
%                 ma_dmaxdev_std(i)   = gen_ma.dmaxdev_std;
%                 edge_dmaxdev_std(i) = gen_edge.dmaxdev_std;
%                 cent_dmaxdev_std(i) = gen_cent.dmaxdev_std;

%                 medaxis_dplus(i)   = KernelDistApproximator(gen_ma.dplus, obs_ma.dplus);
%                 medaxis_dminus(i)  = KernelDistApproximator(gen_ma.dminus, obs_ma.dminus);
%                 edge_dplus(i)      = KernelDistApproximator(gen_edge.dplus, obs_edge.dplus);
%                 edge_dminus(i)     = KernelDistApproximator(gen_edge.dminus, obs_edge.dminus);
%                 centroid_dplus(i)  = KernelDistApproximator(gen_cent.dplus, obs_cent.dplus);
%                 centroid_dminus(i) = KernelDistApproximator(gen_cent.dminus, obs_cent.dminus);

%                 medaxis_result(i)  = get_cdf_result(medaxis_dplus(i), medaxis_dminus(i), data.medaxis_cdf_obs_r(1), data.medaxis_cdf_obs_r(2));
%                 edge_result(i)     = get_cdf_result(edge_dplus(i), edge_dminus(i), data.edge_cdf_obs_r(1), data.edge_cdf_obs_r(2));
%                 centroid_result(i) = get_cdf_result(centroid_dplus(i), centroid_dminus(i), data.centroid_cdf_obs_r(1), data.centroid_cdf_obs_r(2));

                ma_dplus_std(i)    = std(data.medaxis_cdf_gen(1,:));
                ma_dplus_avg(i)    = mean(data.medaxis_cdf_gen(1,:));
                ma_dminus_std(i)   = std(data.medaxis_cdf_gen(2,:));
                ma_dminus_avg(i)   = mean(data.medaxis_cdf_gen(2,:));
                edge_dplus_std(i)  = std(data.edge_cdf_gen(1,:));
                edge_dplus_avg(i)  = mean(data.edge_cdf_gen(1,:));
                edge_dminus_std(i) = std(data.edge_cdf_gen(2,:));
                edge_dminus_avg(i) = mean(data.edge_cdf_gen(2,:));
                cent_dplus_std(i)  = std(data.centroid_cdf_gen(1,:));
                cent_dplus_avg(i)  = mean(data.centroid_cdf_gen(1,:));
                cent_dminus_std(i) = std(data.centroid_cdf_gen(2,:));
                cent_dminus_avg(i) = mean(data.centroid_cdf_gen(2,:));

            end % shape loop
            fprintf('\n');
            results_pval = table(shapes',n,...
                            medaxis_var, medaxis_amd,...
                            centroid_var, centroid_amd,...
                            edge_var, edge_amd, ...
                            medaxis_dplus, medaxis_dminus, medaxis_result, ...
                            centroid_dplus, centroid_dminus, centroid_result, ...
                            edge_dplus, edge_dminus, edge_result);

    %         writetable(results_pval, [out_path patients{p} '_pval_summary_' analysis_conds{c} '.xlsx']);

            results_dstats = table(shapes',n,...
                                    ma_dplus_std, ma_dplus_avg, ma_dminus_std, ma_dminus_avg,...
                                    edge_dplus_std, edge_dplus_avg, edge_dminus_std, edge_dminus_avg,...
                                    cent_dplus_std, cent_dplus_avg, cent_dminus_std, cent_dminus_avg);

    %         writetable(results_dstats, [out_path patients{p} '_dplus_dminus_stats_' analysis_conds{c} '.xlsx']);

            % Store in all subj struct
            all_ma.dist_std(p,:)       = ma_dist_std;
            all_ma.dist_avg(p,:)       = ma_dist_avg;
            all_ma.obs_dist_avg(p,:)   = ma_dist_obs;
            all_ma.obs_dist_var(p,:)   = ma_dist_obs_var;

            all_edge.dist_std(p,:)     = edge_dist_std;
            all_edge.dist_avg(p,:)     = edge_dist_avg;
            all_edge.obs_dist_avg(p,:) = edge_dist_obs;
            all_edge.obs_dist_var(p,:) = edge_dist_obs_var;

            all_cent.dist_std(p,:)     = cent_dist_std;
            all_cent.dist_avg(p,:)     = cent_dist_avg;
            all_cent.obs_dist_avg(p,:) = cent_dist_obs;
            all_cent.obs_dist_var(p,:) = cent_dist_obs_var;

            all_ma.dplus_std(p,:)    = ma_dplus_std;
            all_ma.dplus_avg(p,:)    = ma_dplus_avg;
            all_ma.dminus_std(p,:)   = ma_dminus_std;
            all_ma.dminus_avg(p,:)   = ma_dminus_avg;
            all_ma.obs_dplus(p,:)    = ma_dplus_obs;
            all_ma.obs_dminus(p,:)   = ma_dminus_obs;

            all_edge.dplus_std(p,:)  = edge_dplus_std;
            all_edge.dplus_avg(p,:)  = edge_dplus_avg;
            all_edge.dminus_std(p,:) = edge_dminus_std;
            all_edge.dminus_avg(p,:) = edge_dminus_avg;
            all_edge.obs_dplus(p,:)  = edge_dplus_obs;
            all_edge.obs_dminus(p,:) = edge_dminus_obs;

            all_cent.dplus_std(p,:)  = cent_dplus_std;
            all_cent.dplus_avg(p,:)  = cent_dplus_avg;
            all_cent.dminus_std(p,:) = cent_dminus_std;
            all_cent.dminus_avg(p,:) = cent_dminus_avg;
            all_cent.obs_dplus(p,:)  = cent_dplus_obs;
            all_cent.obs_dminus(p,:) = cent_dminus_obs;

            all_ma.dmaxdev_avg(p,:)   = ma_dmaxdev_avg;
            all_edge.dmaxdev_avg(p,:) = edge_dmaxdev_avg;
            all_cent.dmaxdev_avg(p,:) = cent_dmaxdev_avg;

            all_ma.dmaxdev_std(p,:)   = ma_dmaxdev_std;
            all_edge.dmaxdev_std(p,:) = edge_dmaxdev_std;
            all_cent.dmaxdev_std(p,:) = cent_dmaxdev_std;
        
        end % patient loop 
    
    %     all_column_names = {'Patient',...
    %                         'MedialAxis_AMD_StdDev', 'MedialAxis_AMD_Mean', 'MedialAxis_AMD_Observed',...
    %                         'Centroid_AMD_StdDev', 'Centroid_AMD_Mean', 'Centroid_AMD_Observed',...
    %                         'Edge_AMD_StdDev', 'Edge_AMD_Mean', 'Edge_AMD_Observed',...
    %                         'MedialAxis_Dplus_StdDev', 'MedialAxis_Dplus_Mean', 'MedialAxis_Dplus_Observed',...
    %                         'MedialAxis_Dminus_StdDev', 'MedialAxis_Dminus_Mean', 'MedialAxis_Dminus_Observed',...
    %                         'Centroid_Dplus_StdDev', 'Centroid_Dplus_Mean', 'Centroid_Dplus_Observed',...
    %                         'Centroid_Dminus_StdDev', 'Centroid_Dminus_Mean', 'Centroid_Dminus_Observed',...
    %                         'Edge_Dplus_StdDev', 'Edge_Dplus_Mean', 'Edge_Dplus_Observed',...
    %                         'Edge_Dminus_StdDev', 'Edge_Dminus_Mean', 'Edge_Dminus_Observed',...
    %                         'MedialAxis_DMaxDeviation_Mean', 'Centroid_DMaxDeviation_Mean', 'Edge_DMaxDeviation_Mean',...
    %                         'MedialAxis_DMaxDeviation_StdDev', 'Centroid_DMaxDeviation_StdDev', 'Edge_DMaxDeviation_StdDev'}
    %     
        all_column_names = {'Patient',...
                            'MedialAxis_AMD_Observed','MedialAxis_Var_Observed',...
                            'Centroid_AMD_Observed','Centroid_Var_Observed',...
                            'Edge_AMD_Observed','Edge_Var_Observed'};

    %     results_all = table(patients',...
    %                         all_ma.dist_std, all_ma.dist_avg, all_ma.obs_dist_avg,...
    %                         all_cent.dist_std, all_cent.dist_avg, all_cent.obs_dist_avg,...
    %                         all_edge.dist_std, all_edge.dist_avg, all_edge.obs_dist_avg,...
    %                         all_ma.dplus_std, all_ma.dplus_avg, all_ma.obs_dplus,...
    %                         all_ma.dminus_std, all_ma.dminus_avg, all_ma.obs_dminus,...
    %                         all_cent.dplus_std, all_cent.dplus_avg, all_cent.obs_dplus,...
    %                         all_cent.dminus_std, all_cent.dminus_avg, all_cent.obs_dminus,...
    %                         all_edge.dplus_std, all_edge.dplus_avg, all_edge.obs_dplus,...
    %                         all_edge.dminus_std, all_edge.dminus_avg, all_edge.obs_dminus,...
    %                         all_ma.dmaxdev_avg, all_cent.dmaxdev_avg, all_edge.dmaxdev_avg,...
    %                         all_ma.dmaxdev_std, all_cent.dmaxdev_std, all_edge.dmaxdev_std);

        results_all = table(patients',...
                            all_ma.obs_dist_avg, all_ma.obs_dist_var,...
                            all_cent.obs_dist_avg, all_cent.obs_dist_var,...
                            all_edge.obs_dist_avg, all_edge.obs_dist_var);


        results_all.Properties.VariableNames = all_column_names;
    %   
        if isempty(tasks)
            out_name_all = ['D:\ShapeTapper-Analysis\Stats_allparts_allshapes_' analysis_conds{c} '_amd_var.xlsx'];
        else
            out_name_all = ['D:\ShapeTapper-Analysis\Stats_allparts_allshapes_' analysis_conds{c} '_' tasks{t} '_amd_var.xlsx'];
        end

        writetable(results_all, out_name_all);
    %     
    %     results_allmeans = table(patients',...
    %                         mean(all_ma.dist_std,2), mean(all_ma.dist_avg,2), mean(all_ma.obs_dist_avg,2),...
    %                         mean(all_cent.dist_std,2), mean(all_cent.dist_avg,2), mean(all_cent.obs_dist_avg,2),... 
    %                         mean(all_edge.dist_std,2), mean(all_edge.dist_avg,2), mean(all_edge.obs_dist_avg,2),...                       all_ma.dplus_std, all_ma.dplus_avg,...
    %                         mean(all_ma.dplus_std,2), mean(all_ma.dplus_avg,2), mean(all_ma.obs_dplus,2),...    
    %                         mean(all_ma.dminus_std,2), mean(all_ma.dminus_avg,2), mean(all_ma.obs_dminus,2),...
    %                         mean(all_cent.dplus_std,2), mean(all_cent.dplus_avg,2), mean(all_cent.obs_dplus,2),...
    %                         mean(all_cent.dminus_std,2), mean(all_cent.dminus_avg,2), mean(all_cent.obs_dminus,2),...
    %                         mean(all_edge.dplus_std,2), mean(all_edge.dplus_avg,2), mean(all_edge.obs_dplus,2),...
    %                         mean(all_edge.dminus_std,2), mean(all_edge.dminus_avg,2), mean(all_edge.obs_dminus,2),...
    %                         mean(all_ma.dmaxdev_avg,2), mean(all_cent.dmaxdev_avg,2), mean(all_edge.dmaxdev_avg,2),...    
    %                         mean(all_ma.dmaxdev_std,2), mean(all_cent.dmaxdev_std,2), mean(all_edge.dmaxdev_std,2));

        allmeans_column_names = {'Patient',...
                        'MedialAxis_DPlus_Observed','MedialAxis_DMinus_Observed',...
                        'Centroid_DPlus_Observed','Centroid_DMinus_Observed',...
                        'Edge_DPlus_Observed','Edge_DMinus_Observed'};
    
        results_allmeans = table(patients',...
                            mean(all_ma.obs_dplus,2),mean(all_ma.obs_dminus,2),...
                            mean(all_cent.obs_dplus,2),mean(all_cent.obs_dminus,2),...
                            mean(all_edge.obs_dplus,2),mean(all_edge.obs_dminus,2));
       
        results_allmeans.Properties.VariableNames = allmeans_column_names;

        if isempty(tasks)
            out_name_means = ['D:\ShapeTapper-Analysis\Stats_allparts_shapemeans_' analysis_conds{c} '_dplus_dminus.xlsx'];
        else
            out_name_means = ['D:\ShapeTapper-Analysis\Stats_allparts_shapemeans_' analysis_conds{c} '_' tasks{t} '_dplus_dminus.xlsx'];
        end
    
        writetable(results_allmeans, out_name_means);

    end % task loop

end % condition loop