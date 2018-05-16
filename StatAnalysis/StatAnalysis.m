clear

patients = {'DF','MC',...
            'S01','S02','S03','S04','S06','S07',...
            'S08','S09','S11','S13','S14','S17',...
            'S20','S21','S22','S23','S24'};
% patients = {'DF','MC'};
num_patients = length(patients);

analysis_conds = {'in_shape','bounding_circle','touchpoint_hull'};

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

[all_ma.obs_dist_std, all_ma.obs_dist_avg] = deal(zeros(num_patients,num_shapes));
[all_edge.obs_dist_std, all_edge.obs_dist_avg] = deal(zeros(num_patients,num_shapes));
[all_cent.obs_dist_std, all_cent.obs_dist_avg] = deal(zeros(num_patients,num_shapes));

% Init Dplus/Dminus structs for all subjs
[all_ma.dplus_std, all_ma.dplus_avg, all_ma.dminus_std, all_ma.dminus_avg] = deal(zeros(num_patients,num_shapes));
[all_edge.dplus_std, all_edge.dplus_avg, all_edge.dminus_std, all_edge.dminus_avg] = deal(zeros(num_patients,num_shapes));
[all_cent.dplus_std, all_cent.dplus_avg, all_cent.dminus_std, all_cent.dminus_avg] = deal(zeros(num_patients,num_shapes));

[all_ma.obs_dplus, all_ma.obs_dminus] = deal(zeros(num_patients,num_shapes));
[all_edge.obs_dplus, all_edge.obs_dminus] = deal(zeros(num_patients,num_shapes));
[all_cent.obs_dplus, all_cent.obs_dminus] = deal(zeros(num_patients,num_shapes));


for c=1:length(analysis_conds)
    
    for p=1:length(patients)

        out_path = ['D:\ShapeTapper-Analysis\' patients{p} '\spatial_analysis\'];
        if ~exist(out_path, 'dir')
            mkdir(out_path);
        end
        
        fprintf('Starting %s, %s ', patients{p}, analysis_conds{c});
        
        in_path  = ['D:\ShapeTapper-Analysis\' patients{p} '\spatial_analysis\' analysis_conds{c} '\'];

        n = zeros(num_shapes,1);
        
        [ma_dist_std, ma_dist_avg] = deal(zeros(num_shapes,1));
        [edge_dist_std, edge_dist_avg] = deal(zeros(num_shapes,1));
        [cent_dist_std, cent_dist_avg] = deal(zeros(num_shapes,1));
        
        [ma_dist_obs, edge_dist_obs, cent_dist_obs] = deal(zeros(num_shapes,1));
        
        [medaxis_var, medaxis_amd, edge_var, edge_amd, centroid_var, centroid_amd] = deal(zeros(num_shapes,1));
        [medaxis_dplus, medaxis_dminus, edge_dplus, edge_dminus, centroid_dplus, centroid_dminus] = deal(zeros(num_shapes,1));
        [medaxis_result, edge_result, centroid_result] = deal(cell(num_shapes,1));

        [ma_dplus_std, ma_dplus_avg, ma_dminus_std, ma_dminus_avg] = deal(zeros(num_shapes,1));
        [edge_dplus_std, edge_dplus_avg, edge_dminus_std, edge_dminus_avg] = deal(zeros(num_shapes,1));
        [cent_dplus_std, cent_dplus_avg, cent_dminus_std, cent_dminus_avg] = deal(zeros(num_shapes,1));
        
        [ma_dplus_obs, edge_dplus_obs, cent_dplus_obs] = deal(zeros(num_shapes,1));
        [ma_dminus_obs, edge_dminus_obs, cent_dminus_obs] = deal(zeros(num_shapes,1));
        
        for i=1:num_shapes
            fprintf('.');
            spat_file = [in_path shapes{i} '_Patient_' patients{p} '_spatial_analysis_' analysis_conds{c} '.mat'];
            
            try
                data = load(spat_file);
            catch
                continue
            end

            n(i) = data.n_points;
            
            medaxis_var(i)  = KernelDistApproximator(data.uniform_medaxis_data(1,:)', data.observed_medaxis_data(1)); % variance
            medaxis_amd(i)  = KernelDistApproximator(data.uniform_medaxis_data(3,:)', data.observed_medaxis_data(3)); % amd 
            edge_var(i)     = KernelDistApproximator(data.uniform_edge_data(1,:)', data.observed_edge_data(1)); % variance
            edge_amd(i)     = KernelDistApproximator(data.uniform_edge_data(3,:)', data.observed_edge_data(3)); % amd    
            centroid_var(i) = KernelDistApproximator(data.uniform_centroid_data(1,:)', data.observed_centroid_data(1)); % variance
            centroid_amd(i) = KernelDistApproximator(data.uniform_centroid_data(3,:)', data.observed_centroid_data(3)); % amd    

            ma_dist_std(i)   = std(data.uniform_medaxis_data(3,:),[],2);
            ma_dist_avg(i)   = mean(data.uniform_medaxis_data(3,:),2);
            edge_dist_std(i) = std(data.uniform_edge_data(3,:),[],2); 
            edge_dist_avg(i) = mean(data.uniform_edge_data(3,:),2);
            cent_dist_std(i) = std(data.uniform_centroid_data(3,:),[],2); 
            cent_dist_avg(i) = mean(data.uniform_centroid_data(3,:),2);
            
            ma_dist_obs(i)   = data.observed_medaxis_data(3);
            edge_dist_obs(i) = data.observed_edge_data(3);
            cent_dist_obs(i) = data.observed_centroid_data(3);
            
            [obs_ma, gen_ma]     = parse_cdf_struct(data.medaxis_cdf);
            [obs_edge, gen_edge] = parse_cdf_struct(data.edge_cdf);
            [obs_cent, gen_cent] = parse_cdf_struct(data.centroid_cdf);

            ma_dplus_obs(i)    = obs_ma.dplus;
            ma_dminus_obs(i)   = obs_ma.dminus;
            edge_dplus_obs(i)  = obs_edge.dplus;
            edge_dminus_obs(i) = obs_edge.dminus;
            cent_dplus_obs(i)  = obs_cent.dplus;
            cent_dminus_obs(i) = obs_cent.dminus;
            
            medaxis_dplus(i)   = KernelDistApproximator(gen_ma.dplus, obs_ma.dplus);
            medaxis_dminus(i)  = KernelDistApproximator(gen_ma.dminus, obs_ma.dminus);
            edge_dplus(i)      = KernelDistApproximator(gen_edge.dplus, obs_edge.dplus);
            edge_dminus(i)     = KernelDistApproximator(gen_edge.dminus, obs_edge.dminus);
            centroid_dplus(i)  = KernelDistApproximator(gen_cent.dplus, obs_cent.dplus);
            centroid_dminus(i) = KernelDistApproximator(gen_cent.dminus, obs_cent.dminus);

            medaxis_result(i)  = get_cdf_result(medaxis_dplus(i), medaxis_dminus(i), obs_ma.dplus_r, obs_ma.dminus_r);
            edge_result(i)     = get_cdf_result(edge_dplus(i), edge_dminus(i), obs_edge.dplus_r, obs_edge.dminus_r);
            centroid_result(i) = get_cdf_result(centroid_dplus(i), centroid_dminus(i), obs_cent.dplus_r, obs_cent.dminus_r);

            ma_dplus_std(i)    = gen_ma.dplus_std;
            ma_dplus_avg(i)    = gen_ma.dplus_avg;
            ma_dminus_std(i)   = gen_ma.dminus_std;
            ma_dminus_avg(i)   = gen_ma.dminus_avg;
            edge_dplus_std(i)  = gen_edge.dplus_std;
            edge_dplus_avg(i)  = gen_edge.dplus_avg;
            edge_dminus_std(i) = gen_edge.dminus_std;
            edge_dminus_avg(i) = gen_edge.dminus_avg;
            cent_dplus_std(i)  = gen_cent.dplus_std;
            cent_dplus_avg(i)  = gen_cent.dplus_avg;
            cent_dminus_std(i) = gen_cent.dminus_std;
            cent_dminus_avg(i) = gen_cent.dminus_avg;
            
        end % shape loop
        fprintf('\n');
        results_pval = table(shapes',n,...
                        medaxis_var, medaxis_amd,...
                        centroid_var, centroid_amd,...
                        edge_var, edge_amd, ...
                        medaxis_dplus, medaxis_dminus, medaxis_result, ...
                        centroid_dplus, centroid_dminus, centroid_result, ...
                        edge_dplus, edge_dminus, edge_result);

        writetable(results_pval, [out_path patients{p} '_pval_summary_' analysis_conds{c} '.xlsx']);
        
        results_dstats = table(shapes',n,...
                                ma_dplus_std, ma_dplus_avg, ma_dminus_std, ma_dminus_avg,...
                                edge_dplus_std, edge_dplus_avg, edge_dminus_std, edge_dminus_avg,...
                                cent_dplus_std, cent_dplus_avg, cent_dminus_std, cent_dminus_avg);

        writetable(results_dstats, [out_path patients{p} '_dplus_dminus_stats_' analysis_conds{c} '.xlsx']);

        % Store in all subj struct
        all_ma.dist_std(p,:)       = ma_dist_std;
        all_ma.dist_avg(p,:)       = ma_dist_avg;
        all_ma.obs_dist_avg(p,:)   = ma_dist_obs;
        
        all_edge.dist_std(p,:)     = edge_dist_std;
        all_edge.dist_avg(p,:)     = edge_dist_avg;
        all_edge.obs_dist_avg(p,:) = edge_dist_obs;
        
        all_cent.dist_std(p,:)     = cent_dist_std;
        all_cent.dist_avg(p,:)     = cent_dist_avg;
        all_cent.obs_dist_avg(p,:) = cent_dist_obs;
        
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
        
    end % patient loop
    
    all_column_names = {'Patient',...
                        'MedialAxis_AMD_StdDev', 'MedialAxis_AMD_Mean', 'MedialAxis_AMD_Observed',...
                        'Centroid_AMD_StdDev', 'Centroid_AMD_Mean', 'Centroid_AMD_Observed',...
                        'Edge_AMD_StdDev', 'Edge_AMD_Mean', 'Edge_AMD_Observed',...
                        'MedialAxis_Dplus_StdDev', 'MedialAxis_Dplus_Mean', 'MedialAxis_Dplus_Observed',...
                        'MedialAxis_Dminus_StdDev', 'MedialAxis_Dminus_Mean', 'MedialAxis_Dminus_Observed',...
                        'Centroid_Dplus_StdDev', 'Centroid_Dplus_Mean', 'Centroid_Dplus_Observed',...
                        'Centroid_Dminus_StdDev', 'Centroid_Dminus_Mean', 'Centroid_Dminus_Observed',...
                        'Edge_Dplus_StdDev', 'Edge_Dplus_Mean', 'Edge_Dplus_Observed',...
                        'Edge_Dminus_StdDev', 'Edge_Dminus_Mean', 'Edge_Dminus_Observed' };
    
    results_all = table(patients',...
                        all_ma.dist_std, all_ma.dist_avg, all_ma.obs_dist_avg,...
                        all_cent.dist_std, all_cent.dist_avg, all_cent.obs_dist_avg,...
                        all_edge.dist_std, all_edge.dist_avg, all_edge.obs_dist_avg,...
                        all_ma.dplus_std, all_ma.dplus_avg, all_ma.obs_dplus,...
                        all_ma.dminus_std, all_ma.dminus_avg, all_ma.obs_dminus,...
                        all_cent.dplus_std, all_cent.dplus_avg, all_cent.obs_dplus,...
                        all_cent.dminus_std, all_cent.dminus_avg, all_cent.obs_dminus,...
                        all_edge.dplus_std, all_edge.dplus_avg, all_edge.obs_dplus,...
                        all_edge.dminus_std, all_edge.dminus_avg, all_edge.obs_dminus);
                    
    results_all.Properties.VariableNames = all_column_names;
                    
    writetable(results_all, ['D:\ShapeTapper-Analysis\Stats_allparts_allshapes_' analysis_conds{c} '.xlsx']);
    
    results_allmeans = table(patients',...
                        mean(all_ma.dist_std,2), mean(all_ma.dist_avg,2), mean(all_ma.obs_dist_avg,2),...
                        mean(all_cent.dist_std,2), mean(all_cent.dist_avg,2), mean(all_cent.obs_dist_avg,2),... 
                        mean(all_edge.dist_std,2), mean(all_edge.dist_avg,2), mean(all_edge.obs_dist_avg,2),...                       all_ma.dplus_std, all_ma.dplus_avg,...
                        mean(all_ma.dplus_std,2), mean(all_ma.dplus_avg,2), mean(all_ma.obs_dplus,2),...    
                        mean(all_ma.dminus_std,2), mean(all_ma.dminus_avg,2), mean(all_ma.obs_dminus,2),...
                        mean(all_cent.dplus_std,2), mean(all_cent.dplus_avg,2), mean(all_cent.obs_dplus,2),...
                        mean(all_cent.dminus_std,2), mean(all_cent.dminus_avg,2), mean(all_cent.obs_dminus,2),...
                        mean(all_edge.dplus_std,2), mean(all_edge.dplus_avg,2), mean(all_edge.obs_dplus,2),...
                        mean(all_edge.dminus_std,2), mean(all_edge.dminus_avg,2), mean(all_edge.obs_dminus,2));
    
    results_allmeans.Properties.VariableNames = all_column_names;
                    
    writetable(results_allmeans, ['D:\ShapeTapper-Analysis\Stats_allparts_shapemeans_' analysis_conds{c} '.xlsx']);

end % condition loop