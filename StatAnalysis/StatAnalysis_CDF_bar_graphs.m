clear

patients = {'DF','MC'};
dat_path = 'D:\ShapeTapper-Analysis\';

analysis_conds = {'in_shape','bounding_circle','touchpoint_hull','patient_fitted'};

shapes = {'blake_01','blake_03','blake_04','blake_06','blake_07',...
          'blake_08','blake_09','blake_10','blake_11','blake_12'};

% Groups to compare significance measures
% shape_group_1 = {'blake_04'};%{'blake_06','blake_08'};
% shape_group_2 = {};%{'blake_03','blake_07'};
n_shapes   = length(shapes);
% n_shapes_1 = length(shape_group_1);
% n_shapes_2 = length(shape_group_2);

% Shape numbers only (for graph labels)
% gnums_1 = erase(shape_group_1,"blake_");
% gnums_2 = erase(shape_group_2,"blake_");

ro_fields = {"medaxis_cdf","edge_cdf","centroid_cdf"};
ro_titles = {"Medial Axis","Edge","Centroid"};
ro_ftitles = {"MedialAxis","Edge","Centroid"};
line_types = {'--',':','-.'};

for c=1:length(analysis_conds)
    
    dvals = struct;
    dvals.shape_index = shapes;
    
    for ro=1:length(ro_fields)
        
        out_path = [dat_path '\figures\'];
        
        for p=1:length(patients)
            
            in_path  = [dat_path patients{p} '\spatial_analysis\' analysis_conds{c} '\'];
            
            %%%%%% ALL BLAKE D+/D-s %%%%%%
            dplus_obs   = single(zeros(n_shapes,1));
            dminus_obs  = single(zeros(n_shapes,1));
            dplus_unif  = single(zeros(n_shapes,100000));
            dminus_unif = single(zeros(n_shapes,100000));
            
            for i=1:n_shapes
                dat_file = [in_path shapes{i}  '_Patient_' patients{p} '_spatial_analysis_' analysis_conds{c} '.mat'];
                try
                    data = load(dat_file);
                catch ME
                    warning(strcat("Problem loading ", shapes{i}, " for ", patients{p}, ". Skipping...."));
                    continue
                end
                
                expected_cdf = data.(ro_fields{ro}){1,1}{1,1};
                observed_cdf = data.(ro_fields{ro}){1,1}{1,2};
                uniform_cdf  = data.(ro_fields{ro}){1,1}{1,3};
    %             cdf_axis = data.(ro_fields{ro}){1,2};
    
                % Get D+ and D- values, copy to group averages
                % Observed D+/D-
                dplus_obs(i)  = max(observed_cdf - expected_cdf);
                dminus_obs(i) = min(observed_cdf - expected_cdf);
                % Uniform D+/D-
                dplus_unif(i,:)  = max(uniform_cdf - expected_cdf,[],2)';
                dminus_unif(i,:) = min(uniform_cdf - expected_cdf,[],2)';
            end
            
            dvals.(patients{p}).(ro_fields{ro}).dplus_obs   = dplus_obs;
            dvals.(patients{p}).(ro_fields{ro}).dminus_obs  = dminus_obs;
            dvals.(patients{p}).(ro_fields{ro}).dplus_unif  = dplus_unif;
            dvals.(patients{p}).(ro_fields{ro}).dminus_unif = dminus_unif;
            
        end % patient loop
        
    end % reference object loop
    
    %%%%%%%%%%%%%%%%% MAKE BAR GRAPHS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ro=1:length(ro_fields)
        colors = {'r','b'};
        
        %%%%%%% PLOT D+ BAR GRAPH %%%%%%%%%%%%%%%%
        fig1 = figure; hold on
        ax = gca;
        for p=1:length(patients{p})
            dplus_obs  = dvals.(patients{p}).(ro_fields{ro}).dplus_obs;
            dplus_unif = dvals.(patients{p}).(ro_fields{ro}).dplus_unif;
            dplus_unif_means = mean(dplus_unif,2);
            dplus_unif_stds  = std(dplus_unif,0,2);
            dplus_unif_stes  = std(dplus_unif,0,2)/sqrt(length(dplus_unif));
%             x_axis = linspace(1,length(dplus_unif_means), length(dplus_unif_means));
            x_axis = [1,2,3,4,5,6];
            dplus_obs([1 2 6 7])        = [];
            dplus_unif_means([1 2 6 7]) = [];
            dplus_unif_stds([1 2 6 7])  = [];
            dplus_unif_stes([1 2 6 7])  = [];
%             if strcmp(patients{p}, 'DF')
%                 x_axis(2) = [];
%                 dplus_obs(2)        = [];
%                 dplus_unif_means(2) = [];
%                 dplus_unif_stds(2)  = [];
%                 dplus_unif_stes(2)  = [];
%             end
            e(p) = errorbar(x_axis, dplus_unif_means, dplus_unif_stds,'o');
            e(p).LineWidth = 3;
            e(p).Color = colors{p};
            e(p).Marker = 's';
            e(p).MarkerFaceColor = colors{p};
            s(p)=scatter(x_axis, dplus_obs,180,'filled','+','LineWidth',3);
            s(p).MarkerEdgeColor = colors{p};
            s(p).MarkerFaceColor = colors{p};
        end
        
%         ax.XLim = [0 length(dplus_unif_means)+1];
%         shape_labels = {'','blake_01','blake_03','blake_04','blake_06','blake_07','blake_08','blake_09','blake_10','blake_11','blake_12',''};
        ax.XLim = [0 7];
        shape_labels = {'','blake_04','blake_06','blake_07','blake_10','blake_11','blake_12',''};
        set(gca,'xticklabel',shape_labels,'FontSize', 12);
        ax.TickLabelInterpreter = 'none';
        ax.XAxisLocation = 'bottom';
        ylabel('D+ Value', 'FontSize', 18)
        leg_labels = {'DF - Observed','MC - Observed','DF - Simulated Uniform', 'MC - Simulated Uniform'};
        l = legend([s(1) s(2) e(1), e(2) ], leg_labels, 'FontSize', 15, 'location', 'northwest');
        title(strcat(ro_titles{ro}, " D+ Values, Observed vs Uniform"), 'FontSize', 18);
        set(fig1, 'position',[10 10 1200 800]);
        saveas(fig1,strcat(out_path,shapes{i},"_",analysis_conds{c},"_",ro_ftitles{ro},"_DPlus_BarPlot.png"));
        
        %%%%%%% PLOT D- BAR GRAPH %%%%%%%%%%%%%%%%
        fig2 = figure; hold on
        ax = gca;
        for p=1:length(patients{p})
            dminus_obs  = dvals.(patients{p}).(ro_fields{ro}).dminus_obs;
            dminus_unif = dvals.(patients{p}).(ro_fields{ro}).dminus_unif;
            dminus_unif_means = mean(dminus_unif,2);
            dminus_unif_stds  = std(dminus_unif,0,2);
            dminus_unif_stes  = std(dminus_unif,0,2)/sqrt(length(dminus_unif));
            x_axis = [1,2,3,4,5,6];
            dminus_obs([1 2 6 7])        = [];
            dminus_unif_means([1 2,6 7]) = [];
            dminus_unif_stds([1 2 6 7])  = [];
            dminus_unif_stes([1 2 6 7])  = [];
            e(p) = errorbar(x_axis, dminus_unif_means, dminus_unif_stds,'o');
            e(p).LineWidth = 3;
            e(p).Color = colors{p};
            e(p).Marker = 's';
            e(p).MarkerFaceColor = colors{p};
            s(p)=scatter(x_axis, dminus_obs,150,'filled','+','LineWidth',2);
            s(p).MarkerEdgeColor = colors{p};
            s(p).MarkerFaceColor = colors{p};
        end
        
%         ax.XLim = [0 length(dplus_unif_means)+1];
%         shape_labels = {'','blake_01','blake_03','blake_04','blake_06','blake_07','blake_08','blake_09','blake_10','blake_11','blake_12',''};
        ax.XLim = [0 7];
        shape_labels = {'','blake_04','blake_06','blake_07','blake_10','blake_11','blake_12',''};
        set(gca,'xticklabel',shape_labels,'FontSize', 12);
        ax.TickLabelInterpreter = 'none';
        ax.XAxisLocation = 'bottom';
        ylabel('D+ Value', 'FontSize', 18)
        leg_labels = {'DF - Observed','MC - Observed','DF - Simulated Uniform', 'MC - Simulated Uniform'};
        l = legend([s(1) s(2) e(1), e(2) ], leg_labels, 'FontSize', 15, 'location', 'southwest');
        title(strcat(ro_titles{ro}, " D- Values, Observed vs Uniform"), 'FontSize', 18);
        set(fig2, 'position',[10 10 1200 800]);
        saveas(fig2,strcat(out_path,shapes{i},"_",analysis_conds{c},"_",ro_ftitles{ro},"_DMinus_BarPlot.png"));
    end % RO loop 2 (plots)
    
end % condition loop