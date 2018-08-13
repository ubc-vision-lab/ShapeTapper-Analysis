clear

% patients = {'MC1','MC2'};
patients = {'DF'};
        
img_names = ["blake_01","blake_03","blake_04","blake_06","blake_07",...
             "blake_08","blake_09","blake_10","blake_11","blake_12"];
        
% tasks = {'Simultaneous_2AFC','Sequential_2AFC','Oddball','N_Back',''};
% task_keys = {'lr','nb','ob','solo'};

in_path  = 'C:\Users\Visionlab\Development\ShapeTapper-Analysis\';
out_path = 'D:\ShapeTapper-Analysis\';


for p=1:length(patients)
    
    dat_path = [in_path patients{p} '\aggregated_observations\raw_touchpoint_data\'];
    files = dir(dat_path);
    fnames = {files.name};
    
    fileID = fopen();
    C1 = textscan(fileID,'%s %s %f %f %s','headerLines',1);
    fclose(fileID);

%     dat = table(C1{1}, C1{2}, [cell2mat(C1(3)) cell2mat(C1(4))], 'VariableNames',{'OriginFileName','ImageName','XY'});

    fnames = C1{1};  % odd matrix
    inames = C1{5};  % even matrix
    xs = C1{2}; 
    ys = num2cell(C1{3}); 
    
    dat = table(fnames, inames, [xs ys], 'VariableNames',{'OriginFileName','ImageName','XY'});
    
    for t=1:length(tasks)
        if strcmp(tasks{t}, '')
            task = true(length(fnames),1);
        else
            task = contains(fnames, task_keys{t});
        end
        
        dat = table(fnames(task), inames(task), [xs(task) ys(task)], 'VariableNames',{'OriginFileName','ImageName','XY'});
  
        % fileID = fopen("..\MC2\aggregated_observations\Overall_Results.txt");
        % C2 = textscan(fileID,'%s %s %f %f','headerLines',1);
        % fclose(fileID);
        %
        % dat2 = table(C2{1}, C2{2}, [cell2mat(C2(3)) cell2mat(C2(4))], 'VariableNames',{'OriginFileName','ImageName','XY'});
        % dat2 = dat2(~strcmp(dat2.OriginFileName, 's_solo_blake_10b36p_3m_Y8hZ'),:);
        % 
        % dat = [dat; dat2];

        for i=1:length(img_names)

            shape_dat = dat(strcmp(dat.ImageName, img_names(i)),:);
            img_dataset = shape_dat.XY;

            out_dir = [out_path patients{p} '\observed_touchpoints\'];

            if ~exist(out_dir, 'dir')
                mkdir(out_dir);
            end
            
            if strcmp(tasks{t}, '')
                output_name = [out_dir img_names{i} '_Patient_' patients{p} '_observed_touchpoints.mat'];
            else
                output_name = [out_dir img_names{i} '_Patient_' patients{p} '_observed_touchpoints_' tasks{t} '.mat'];
            end
            
            save(output_name,'img_dataset');

        end %shape loop
    end % task loop
end %patient loop