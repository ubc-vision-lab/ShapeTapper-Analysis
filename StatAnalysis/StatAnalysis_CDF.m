clear

patient = 'MC';
dat_path = 'D:\ShapeTapper-Analysis\';

analysis_conds = {'in_shape'};%{'in_shape','bounding_circle','touchpoint_hull','patient_fitted'};

shapes = {'blake_01','blake_03','blake_04','blake_06','blake_07','blake_08','blake_09','blake_10','blake_11','blake_12'};

% Groups to compare significance measures
shape_group_1 = {'blake_04'};%{'blake_06','blake_08'};
shape_group_2 = {};%{'blake_03','blake_07'};
n_shapes   = length(shapes);
n_shapes_1 = length(shape_group_1);
n_shapes_2 = length(shape_group_2);

% Shape numbers only (for graph labels)
gnums_1 = erase(shape_group_1,"blake_");
gnums_2 = erase(shape_group_2,"blake_");

ro_fields = {"medaxis_cdf"};
line_types = {'--',':','-.'};
      
for c=1:length(analysis_conds)
    
    pvals = struct;
    
    in_path  = [dat_path patient '\distance_analysis\' analysis_conds{c} '\'];

    for ro=1:length(ro_fields)
        
        %%%%%% ALL BLAKE CDFs %%%%%%
        all_exp_cdfs = figure; hold on
        for i=1:n_shapes
            dat_file = [in_path shapes{i} '_analysis.mat'];
            data = load(dat_file);

            expected_cdf = data.(ro_fields{ro}){1,1}{1,1};
            observed_cdf = data.(ro_fields{ro}){1,1}{1,2};
%             cdf_axis = data.(ro_fields{ro}){1,2};

            plot(expected_cdf, line_types{mod(i,3)+1}, 'LineWidth', 2);
        end
        t = title(strcat("Expected CDF (", erase(ro_fields{ro},"_cdf"),") All Shapes"));
        set(t,'Interpreter', 'none')
        leg=legend(shapes);
        set(leg,'Interpreter', 'none')
        saveas(all_exp_cdfs, strjoin(["figures\" string(patient), ro_fields{ro}, "all_expected.png"], "_"));

        %%%%%% GROUP BLAKE CDFs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%% GROUP 1 %%%%%%
        exp_cdf_g1  = single(zeros(n_shapes_1,1000));
        obs_cdf_g1  = single(zeros(n_shapes_1,1000));
        dplus_o_g1  = single(zeros(n_shapes_1,2));
        dminus_o_g1 = single(zeros(n_shapes_1,2));
        dplus_u_g1  = single(zeros(n_shapes_1,100000));
        dminus_u_g1 = single(zeros(n_shapes_1,100000));
        
        group_avg_cdf = figure; hold on
        for i=1:n_shapes_1
            dat_file = [in_path shape_group_1{i} '_analysis.mat'];
            data = load(dat_file);
            
%             cdf_axis = data.(ro_fields{ro}){1,2};
            expected_cdf = data.(ro_fields{ro}){1,1}{1,1};
            observed_cdf = data.(ro_fields{ro}){1,1}{1,2};
            exp_cdf_g1(i,:) = expected_cdf;
            obs_cdf_g1(i,:) = observed_cdf;
            uniform_cdf = data.(ro_fields{ro}){1,1}{1,3};
            
            % Get D+ and D- values, copy to group averages
            % Observed D+/D-
            [dplus_obs,dplus_obs_r]   = max(observed_cdf - expected_cdf);
            [dminus_obs,dminus_obs_r] = min(observed_cdf - expected_cdf);
            dplus_o_g1(i,:)  = [dplus_obs,dplus_obs_r];
            dminus_o_g1(i,:) = [dminus_obs,dminus_obs_r];
            % Uniform D+/D-
            dplus_unif  = max(uniform_cdf - expected_cdf,[],2);
            dminus_unif = min(uniform_cdf - expected_cdf,[],2);
            dplus_u_g1(i,:)  = dplus_unif;
            dminus_u_g1(i,:) = dminus_unif;
            
            plot(expected_cdf, [line_types{i} 'r']);
        end
        plot(mean(exp_cdf_g1,1), '-r', 'LineWidth',2);
        
        %%%%%% GROUP 2 %%%%%%
        exp_cdf_g2 = single(zeros(n_shapes_2,1000));
        obs_cdf_g2 = single(zeros(n_shapes_2,1000));
        dplus_o_g2  = single(zeros(n_shapes_2,2));
        dminus_o_g2 = single(zeros(n_shapes_2,2));
        dplus_u_g2  = single(zeros(n_shapes_2,100000));
        dminus_u_g2 = single(zeros(n_shapes_2,100000));
        
        for i=1:n_shapes_2
            dat_file = [in_path shape_group_2{i} '_analysis.mat'];
            data = load(dat_file);
            
%             cdf_axis = data.(ro_fields{ro}){1,2};
            expected_cdf = data.(ro_fields{ro}){1,1}{1,1};
            observed_cdf = data.(ro_fields{ro}){1,1}{1,2};
            exp_cdf_g2(i,:) = expected_cdf;
            obs_cdf_g2(i,:) = observed_cdf;
            uniform_cdf = data.(ro_fields{ro}){1,1}{1,3};
            
            % Get D+ and D- values, copy to group averages
            % Observed D+/D-
            [dplus_obs,dplus_obs_r]   = max(observed_cdf - expected_cdf);
            [dminus_obs,dminus_obs_r] = min(observed_cdf - expected_cdf);
            dplus_o_g2(i,:)  = [dplus_obs,dplus_obs_r];
            dminus_o_g2(i,:) = [dminus_obs,dminus_obs_r];
            % Uniform D+/D-
            dplus_unif  = max(uniform_cdf - expected_cdf,[],2);
            dminus_unif = min(uniform_cdf - expected_cdf,[],2);
            dplus_u_g2(i,:)  = dplus_unif;
            dminus_u_g2(i,:) = dminus_unif;

            plot(expected_cdf, [line_types{i} 'b']);
        end
        plot(mean(exp_cdf_g2,1), '-b', 'LineWidth',2);
        title(strcat("Expected CDF (", erase(ro_fields{ro},"_cdf"),") Averages by Group"));
        leg_txt = [shape_group_1 {strcat("Mean 1")} shape_group_2 {strcat("Mean 2")}];
        leg=legend(leg_txt);
        set(leg,'Interpreter', 'none')
        saveas(group_avg_cdf, strjoin(["figures\" string(patient), ro_fields{ro}, "group_averages.png"], "_"));
        
        %%%%% PLOT GROUP 1 & 2 EXPECTED VS OBSERVED %%%%%%
        exp_obs_cdfs = figure; hold on
        t = title(strcat("Expected vs Observed CDF (", erase(ro_fields{ro},"_cdf"),") by Group"));
        set(t,'Interpreter', 'none')
        plot(mean(exp_cdf_g1,1), '-r', 'LineWidth',2);
        plot(mean(obs_cdf_g1,1), '--c', 'LineWidth',2);
        plot(mean(exp_cdf_g2,1), '-b', 'LineWidth',2);
        plot(mean(obs_cdf_g2,1), '--m', 'LineWidth',2);
        leg=legend("Group 1 Expected","Group 1 Observed",...
                    "Group 2 Expected","Group 2 Observed"); 
        set(leg,'Interpreter', 'none')
        ht = text(0.95*size(exp_cdf_g1,2), 0.25, ...
                  {'{\color{red} o } Group 1: 04, 06, 08', ...
                   '{\color{blue} o } Group 2: 03, 07, 10'},...
                  'HorizontalAlignment' ,'right','EdgeColor', 'k');
        saveas(exp_obs_cdfs, strjoin(["figures\"  string(patient), ro_fields{ro}, "expected_vs_observed.png"], "_"));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        %%%%% Test significance across groups %%%%%%
        pvals.dplus.(ro_fields{ro}) = zeros(2,2);
        dplus_1to1 = KernelDistApproximator(mean(dplus_u_g1)', mean(dplus_o_g1(:,1)));
        dplus_2to2 = KernelDistApproximator(mean(dplus_u_g2)', mean(dplus_o_g2(:,1)));
        dplus_1to2 = KernelDistApproximator(mean(dplus_u_g2)', mean(dplus_o_g1(:,1)));
        dplus_2to1 = KernelDistApproximator(mean(dplus_u_g1)', mean(dplus_o_g2(:,1)));
        pvals.dplus.(ro_fields{ro})(1,1) = dplus_1to1;
        pvals.dplus.(ro_fields{ro})(2,2) = dplus_2to2;
        pvals.dplus.(ro_fields{ro})(1,2) = dplus_1to2;
        pvals.dplus.(ro_fields{ro})(2,1) = dplus_2to1;
        
        pvals.dminus.(ro_fields{ro}) = zeros(2,2);
        dminus_1to1 = KernelDistApproximator(mean(dminus_u_g1)', mean(dminus_o_g1(:,1)));
        dminus_2to2 = KernelDistApproximator(mean(dminus_u_g2)', mean(dminus_o_g2(:,1)));
        dminus_1to2 = KernelDistApproximator(mean(dminus_u_g2)', mean(dminus_o_g1(:,1)));
        dminus_2to1 = KernelDistApproximator(mean(dminus_u_g1)', mean(dminus_o_g2(:,1)));
        pvals.dminus.(ro_fields{ro})(1,1) = dminus_1to1;
        pvals.dminus.(ro_fields{ro})(2,2) = dminus_2to2;
        pvals.dminus.(ro_fields{ro})(1,2) = dminus_1to2;
        pvals.dminus.(ro_fields{ro})(2,1) = dminus_2to1;
        
        pvals.results.(ro_fields{ro}) = cell(2,2);
        pvals.results.(ro_fields{ro})(1,1) = get_cdf_result(dplus_1to1, dminus_1to1, sum(dplus_o_g1(:,2)), sum(dminus_o_g1(:,2)));
        pvals.results.(ro_fields{ro})(2,2) = get_cdf_result(dplus_2to2, dminus_2to2, sum(dplus_o_g2(:,2)), sum(dminus_o_g2(:,2)));
        pvals.results.(ro_fields{ro})(1,2) = get_cdf_result(dplus_1to2, dminus_1to2, sum(dplus_o_g1(:,2)), sum(dminus_o_g1(:,2)));
        pvals.results.(ro_fields{ro})(2,1) = get_cdf_result(dplus_2to1, dminus_2to1, sum(dplus_o_g2(:,2)), sum(dminus_o_g2(:,2)));
        
    end % end reference object loop
 
end % end condition loop