clear;close all;

patients = {'S01','S02','S03','S04','S06','S07',...
            'S08','S09','S11','S13','S14','S17',...
            'S20','S21','S22','S23','S24',...
            'S25','S26','S27',...
            'S28','S29','S30',...
            'S31','S32','S33',...
            'DF','MC'};
% patients = {'DF','MC'};
num_patients = length(patients);

shapes = {'blake_01','blake_03','blake_04','blake_06','blake_07',...
          'blake_08','blake_09','blake_10','blake_11','blake_12'};
num_shapes = length(shapes);

%%%%%%%%%%%%%%%% PLOT BY PATIENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for p=1:length(patients)

    out_path = ['D:\ShapeTapper-Analysis\bayesian_analysis\'];
    if ~exist(out_path, 'dir')
        mkdir(out_path);
    end
    
    bayes_fig = figure; hold on
    t = title([patients{p} ' - All Shapes']);
    set(t,'Interpreter', 'none')
    
    for i=1:num_shapes
         
        in_path  = ['D:\ShapeTapper-Analysis\' patients{p} '\bayesian_analysis\'];
        
        bayes_file = [in_path shapes{i} '_Patient_' patients{p} '_bayes.mat'];

        try
            data = load(bayes_file);
        catch
            continue
        end

        mdax_prob = data.medaxis_p / data.null_p;
        mdax_solo_prob = data.medaxis_solo_p / data.null_p;
        edge_prob = data.edge_p / data.null_p;
        cent_prob = data.centroid_p / data.null_p;
        
        xaxis_mm = (1:60)*0.2482;
        
        plot(xaxis_mm, ones(1, 60),'k', 'LineWidth', 2)
        plot(xaxis_mm, edge_prob(1:60), 'r')
        plot(xaxis_mm, mdax_solo_prob(1:60), 'b');
%         plot(cent_prob(1:60), 'g')
    end
    f=get(gca,'Children');
    
    legend(f([1:3]),{'Medial Axis - Centroid',' All Edge','Null Hypothesis'});

    xticks(0:15)
    ylim([0 8])
    grid on
    set(gca, 'gridlinestyle', '--')
    xlabel("StdDev of Gauss smoothing in mm")
    ylabel("Bayes Factor (1 = NULL)")

    saveas(bayes_fig, strjoin([out_path, patients{p}, "_bayesian_probs.png"], ""));
end



%%%%%%%%%%%%%%%% PLOT BY SHAPE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:num_shapes

    out_path = ['D:\ShapeTapper-Analysis\bayesian_analysis\'];
    if ~exist(out_path, 'dir')
        mkdir(out_path);
    end
    
    bayes_fig = figure; hold on
    t = title([shapes{i} ' - All Patients']);
    set(t,'Interpreter', 'none')
    
    for p=1:length(patients)
         
        in_path  = ['D:\ShapeTapper-Analysis\' patients{p} '\bayesian_analysis\'];
        
        bayes_file = [in_path shapes{i} '_Patient_' patients{p} '_bayes.mat'];

        try
            data = load(bayes_file);
        catch
            continue
        end

        mdax_prob = data.medaxis_p / data.null_p;
        mdax_solo_prob = data.medaxis_solo_p / data.null_p;
        edge_prob = data.edge_p / data.null_p;
        cent_prob = data.centroid_p / data.null_p;
        
        xaxis_mm = (1:60)*0.2482;
        
        plot(xaxis_mm, ones(1, 60),'k', 'LineWidth', 2)
        plot(xaxis_mm, edge_prob(1:60), 'r')
        
        if strcmp(patients{p}, 'MC')
            plot(xaxis_mm, mdax_solo_prob(1:60), 'g', 'LineWidth', 2);
        elseif strcmp(patients{p}, 'DF')
            plot(xaxis_mm, mdax_solo_prob(1:60), 'm', 'LineWidth', 2);
        else
            plot(xaxis_mm, mdax_solo_prob(1:60), 'b');
        end
        
%         plot(cent_prob(1:60), 'g')
    end
    f=get(gca,'Children');

    if i == 2 % blake_03 has no data for DF
        legend(f([1 7:9]),{'MC (Medial Axis - Centroid)','Control (Medial Axis - Centroid)','All Edge','Null Hypothesis'});
    else
        legend(f([1 4 7:9]),{'DF (Medial Axis - Centroid)', 'MC (Medial Axis - Centroid)','Control (Medial Axis - Centroid)','All Edge','Null Hypothesis'});
    end

    xticks(0:15)
    ylim([0 8])
    grid on
    set(gca, 'gridlinestyle', '--')
    xlabel("StdDev of Gauss smoothing in mm")
    ylabel("Bayes Factor (1 = NULL)")

    saveas(bayes_fig, strjoin([out_path, patients{p}, "_bayesian_probs.png"], ""));
end