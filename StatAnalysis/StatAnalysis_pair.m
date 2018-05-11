clear

% patients = {'DF'};
patients = {'MC'};

analysis_conds = {'in_shape','bounding_circle','touchpoint_hull',...
                  'patient_fitted','normal_distribution'};

% DF
% shapes = {'blake_04_to_blake_07','blake_07_to_blake_04',...
%           'blake_11_to_blake_07','blake_07_to_blake_11',...
%           'blake_10_to_blake_07','blake_07_to_blake_10',...
%           'blake_04_to_blake_10','blake_10_to_blake_04'};

% MC      
shapes = {'blake_03_to_blake_07','blake_07_to_blake_03',...
          'blake_04_to_blake_07','blake_07_to_blake_04',...
          'blake_07_to_blake_10','blake_10_to_blake_07',...
          'blake_04_to_blake_10','blake_10_to_blake_04'}; 
      
shape_labels = strrep(shapes,'_',' ');

num_shapes = length(shapes);

for p=1:length(patients)
    
    out_path = ['D:\ShapeTapper-Analysis\' patients{p} '\spatial_analysis\'];
    
    for c=1:length(analysis_conds)
        
        fprintf('Starting %s, %s ', patients{p}, analysis_conds{c});
        
        in_path  = ['D:\ShapeTapper-Analysis\' patients{p} '\spatial_analysis\' analysis_conds{c} '\'];

        n = zeros(num_shapes,1);
        [medaxis_var, medaxis_amd, edge_var, edge_amd, centroid_var, centroid_amd] = deal(zeros(num_shapes,1));
        [medaxis_dplus, medaxis_dminus, edge_dplus, edge_dminus, centroid_dplus, centroid_dminus] = deal(zeros(num_shapes,1));
        [medaxis_result, edge_result, centroid_result] = deal(cell(num_shapes,1));

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

            [obs_ma_dplus, obs_ma_dplus_r, obs_ma_dminus, obs_ma_dminus_r, gen_ma_dplus, gen_ma_dminus] = parse_cdf_struct(data.medaxis_cdf);
            [obs_edge_dplus, obs_edge_dplus_r, obs_edge_dminus, obs_edge_dminus_r, gen_edge_dplus, gen_edge_dminus] = parse_cdf_struct(data.edge_cdf);
            [obs_cent_dplus, obs_cent_dplus_r, obs_cent_dminus, obs_cent_dminus_r, gen_cent_dplus, gen_cent_dminus] = parse_cdf_struct(data.centroid_cdf);

            medaxis_dplus(i)   = KernelDistApproximator(gen_ma_dplus, obs_ma_dplus);
            medaxis_dminus(i)  = KernelDistApproximator(gen_ma_dminus, obs_ma_dminus);
            edge_dplus(i)      = KernelDistApproximator(gen_edge_dplus, obs_edge_dplus);
            edge_dminus(i)     = KernelDistApproximator(gen_edge_dminus, obs_edge_dminus);
            centroid_dplus(i)  = KernelDistApproximator(gen_cent_dplus, obs_cent_dplus);
            centroid_dminus(i) = KernelDistApproximator(gen_cent_dminus, obs_cent_dminus);

            medaxis_result(i)  = get_cdf_result(medaxis_dplus(i), medaxis_dminus(i), obs_ma_dplus_r, obs_ma_dminus_r);
            edge_result(i)     = get_cdf_result(edge_dplus(i), edge_dminus(i), obs_edge_dplus_r, obs_edge_dminus_r);
            centroid_result(i) = get_cdf_result(centroid_dplus(i), centroid_dminus(i), obs_cent_dplus_r, obs_cent_dminus_r);

        end % shape loop
        fprintf('\n');
        results = table(shape_labels',n,...
                        medaxis_var, medaxis_amd,...
                        centroid_var, centroid_amd,...
                        edge_var, edge_amd, ...
                        medaxis_dplus, medaxis_dminus, medaxis_result, ...
                        centroid_dplus, centroid_dminus, centroid_result, ...
                        edge_dplus, edge_dminus, edge_result);

        writetable(results, [out_path patients{p} '_stat_summary_' analysis_conds{c} '.xlsx']);

    end % condition loop
end % patients{p} loop