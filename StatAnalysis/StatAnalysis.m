patient = 'MC';

analysis_conds = {'in_shape','bounding_circle','touchpoint_hull'};

shapes = {'blake_01','blake_04','blake_06','blake_07','blake_08','blake_09','blake_10','blake_11','blake_12',...
          'solo3','solo5','solo6','solo7','solo9','solo10','solo11','solo12'};

for c=1:3

    in_path  = ['D:\ShapeTapper-Analysis\' patient '\distance_analysis\' analysis_conds{c} '/'];
    out_path = ['D:\ShapeTapper-Analysis\' patient '\distance_analysis\'];
    
    num_shapes = length(shapes);
    
    [medaxis_var, medaxis_amd, edge_var, edge_amd, centroid_var, centroid_amd] = deal(zeros(num_shapes,1));
    [medaxis_dplus, medaxis_dminus, edge_dplus, edge_dminus, centroid_dplus, centroid_dminus] = deal(zeros(num_shapes,1));
    [medaxis_result, edge_result, centroid_result] = deal(cell(num_shapes,1));
    
    for i=1:num_shapes
        
        spat_file = [in_path shapes{i} '_analysis.mat'];
        data = load(spat_file);

        if isfield(data,'observed_medaxis_var')
            medaxis_var(i) = KernelDistApproximator(data.generated_medaxis_data(1,:)', data.observed_medaxis_var(1)); % variance
            medaxis_amd(i) = KernelDistApproximator(data.generated_medaxis_data(3,:)', data.observed_medaxis_var(3)); % amd
        elseif isfield(data,'observed_medaxis_data')
            medaxis_var(i) = KernelDistApproximator(data.generated_medaxis_data(1,:)', data.observed_medaxis_data(1)); % variance
            medaxis_amd(i) = KernelDistApproximator(data.generated_medaxis_data(3,:)', data.observed_medaxis_data(3)); % amd
        end
           
        edge_var(i) = KernelDistApproximator(data.generated_edge_data(1,:)', data.observed_edge_data(1)); % variance
        edge_amd(i) = KernelDistApproximator(data.generated_edge_data(3,:)', data.observed_edge_data(3)); % amd    
        centroid_var(i) = KernelDistApproximator(data.generated_centroid_data(1,:)', data.observed_centroid_data(1)); % variance
        centroid_amd(i) = KernelDistApproximator(data.generated_centroid_data(3,:)', data.observed_centroid_data(3)); % amd    

        gen_ma_cdf_unif = data.medaxis_cdf{1,1}{1,1};
        gen_ma_cdf_obs  = data.medaxis_cdf{1,1}{1,2};
        obs_ma_dplus = gen_ma_cdf_obs - gen_ma_cdf_unif;
        gen_ma_dplus = data.medaxis_cdf{1,1}
        gen_ma_dminus
        medaxis_dplus(i) = KernelDistApproximator(data.medaxis_cdf.{1,1}generated_medaxis_dplus', data.observed_medaxis_dplus(1));
        medaxis_dminus(i) = KernelDistApproximator(data.generated_medaxis_dminus', data.observed_medaxis_dminus(1));
        edge_dplus(i) = KernelDistApproximator(data.generated_edge_dplus', data.observed_edge_dplus(1));
        edge_dminus(i) = KernelDistApproximator(data.generated_edge_dminus', data.observed_edge_dminus(1));
        centroid_dplus(i) = KernelDistApproximator(data.generated_centroid_dplus', data.observed_centroid_dplus(1));
        centroid_dminus(i) = KernelDistApproximator(data.generated_centroid_dminus', data.observed_centroid_dminus(1));

        if medaxis_dplus(i) > 0.95
            if medaxis_dminus(i) >= 0.05
                medaxis_result(i) = {'N'};
            elseif medaxis_dminus(i) < 0.05 
                if data.observed_medaxis_dplus(2) <= data.observed_medaxis_dminus(2)
                    medaxis_result(i) = {'NF'};
                else 
                    medaxis_result(i) = {'M'};
                end
            end
        elseif medaxis_dminus(i) < 0.05
            medaxis_result(i) = {'F'};
        else 
            medaxis_result(i) = {'Ø'};
        end
        
        if edge_dplus(i) > 0.95
            if edge_dminus(i) >= 0.05
                edge_result(i) = {'N'};
            elseif edge_dminus(i) < 0.05 
                if data.observed_edge_dplus(2) <= data.observed_edge_dminus(2)
                    edge_result(i) = {'NF'};
                else 
                    edge_result(i) = {'M'};
                end
            end
        elseif edge_dminus(i) < 0.05
            edge_result(i) = {'F'};
        else 
            edge_result(i) = {'Ø'};
        end
        
        if centroid_dplus(i) > 0.95
            if centroid_dminus(i) >= 0.05
                centroid_result(i) = {'N'};
            elseif centroid_dminus(i) < 0.05 
                if data.observed_centroid_dplus(2) <= data.observed_centroid_dminus(2)
                    centroid_result(i) = {'NF'};
                else 
                    centroid_result(i) = {'M'};
                end
            end
        elseif centroid_dminus(i) < 0.05
            centroid_result(i) = {'F'};
        else 
            centroid_result(i) = {'Ø'};
        end
        
%         [h,p,stats] = kstest2(data.generated_medaxis_dplus,data.generated_centroid_dplus);
%         [p,h,stats] = ranksum(data.generated_centroid_dplus,data.generated_medaxis_dplus);
%         
%         fprintf('%s : %d %f %d\n', shapes{i}, h, p, stats.ranksum);
%         
    end % end shape loop
    
    results = table(shapes',medaxis_var, medaxis_amd, centroid_var, centroid_amd, edge_var, edge_amd, ...
                    medaxis_dplus, medaxis_dminus, medaxis_result, ...
                    centroid_dplus, centroid_dminus, centroid_result, ...
                    edge_dplus, edge_dminus, edge_result);
                
%     save(['stat_summary_' analysis_conds{c} '.mat'], 'results');
    writetable(results, [out_path 'stat_summary_' analysis_conds{c} '.xlsx']);
%     
end % end condition loop