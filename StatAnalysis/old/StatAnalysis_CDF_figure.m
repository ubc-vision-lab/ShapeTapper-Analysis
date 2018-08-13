clear

patients = {'DF','MC'};

dat_path = 'D:\ShapeTapper-Analysis\';

analysis_conds = {'in_shape'};%{'in_shape','bounding_circle','touchpoint_hull','patient_fitted'};

shapes = {'blake_01','blake_03','blake_04','blake_06','blake_07','blake_08','blake_09','blake_10','blake_11','blake_12'};
n_shapes   = length(shapes);

ro_fields = {"medaxis_cdf"};
line_types = {'--',':','-.'};

for p=1:length(patients)
    
    out_path = [dat_path patients{p} '\figures\' 'cdf_figs\'];
    
    for c=1:length(analysis_conds)

        in_path  = [dat_path patients{p} '\spatial_analysis\' analysis_conds{c} '\'];

        for ro=1:length(ro_fields)

            for i=1:n_shapes
                dat_file = [in_path shapes{i} '_Patient_' patients{p} '_spatial_analysis_' analysis_conds{c} '.mat'];

                try
                    data = load(dat_file);
                catch
                    continue
                end
                
                expected_cdf = data.(ro_fields{ro}){1,1}{1,1};
                observed_cdf = data.(ro_fields{ro}){1,1}{1,2};
                uniform_cdf = data.(ro_fields{ro}){1,1}{1,3};

                % Get D+ and D- values, copy to group averages
                % Observed D+/D-
                [dplus_obs,dplus_obs_r]   = max(observed_cdf - expected_cdf);
                [dminus_obs,dminus_obs_r] = min(observed_cdf - expected_cdf);
                % Uniform D+/D-
                dplus_unif  = max(uniform_cdf - expected_cdf,[],2);
                dminus_unif = min(uniform_cdf - expected_cdf,[],2);

                exp_obs_cdfs = figure; hold on
                set(exp_obs_cdfs, 'Visible', 'off');
                t = title(strcat("Expected vs Observed CDF (", shapes{i}, ")"));
                set(t,'Interpreter', 'none')
                idx_one = 50*(ceil(find(expected_cdf==1)/50.));
                x_ax = linspace(0,100,idx_one(1));
                plot(x_ax,expected_cdf(1:idx_one(1))*100, '-b', 'LineWidth',2);
                plot(x_ax,observed_cdf(1:idx_one(1))*100, '--c', 'LineWidth',2);
                plot([33 33], [0 100], '--r', 'LineWidth',1);
                plot([67 67], [0 100], '--r', 'LineWidth',1);
                xtickformat('percentage'); ytickformat('percentage');
                xlabel('Percentage of total shape canvas'); % y-axis label
                ylabel('Percentage of touchpoints in region'); % y-axis label
                leg=legend("Expected Value (shape area)","Observed Value (touchpoints)","Region=33%","Region=67%");
                set(leg,'Interpreter', 'none', 'Location','east')
                saveas(exp_obs_cdfs, strcat(out_path, shapes{i},"_Patient_",patients{p},"_cdf.png"));

            end % END SHAPE LOOP
        end % end reference object loop
    end % end condition loop
end