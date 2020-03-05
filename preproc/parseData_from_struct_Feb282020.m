clear

patients = {'S01','S02','S03','S04','S05','S06',...
            'S07','S08','S09','S10','S11','S12',...
            'S13','S14','S15','S16','S17','S18',...
            'S19','S20','S21','S22','S23','S24',...
            'S25','S26','S27','S28','S29','S30',...
            'S31','S32','S33','S34','S35','S36',...
            'S37','S38','S39','S40','S41','S42'};
% patients = {'S25','S26','S27','S28','S29','S30',...
%             'S31','S32','S33','S34','S35','S36',... 
%             'S37','S38','S39','S40','S41','S42',...
%             'S43','S44','S45','S46','S47','S48',...
%             'S49','S50','S51','S52','S53','S54'};
        
% patient_nums = [1:54];


img_names = ["blake_01","blake_04","blake_06","blake_07",...
             "blake_08","blake_10","blake_11","blake_12"];
       
tasks = {'1','2'};
% tasks = {'Simultaneous_2AFC','Sequential_2AFC','Oddball','N_Back'};
% task_keys = {'lr','nb','ob','solo'};

% in_path  = 'D:\ShapeTapper-Analysis\';
out_path = 'C:\ShapeTapper-Analysis\';

% load('tps.mat') % add whichever MAT file contains "A" here.

% % Get participant ID list, sort alphanumerically ignoring case
% PIDs = unique(A(2:end,1));
% [~,idx]=sort(upper(PIDs));
% PIDs = PIDs(idx);
% save('PIDs.mat',PIDs)

% load('PIDs.mat')

for p=1:length(patients)
    
%     fileID = fopen([in_path patients{p} '\' num2str(patient_nums(p)) '_Overall_Results.txt']);
%     C1 = textscan(fileID,'%s %f %f %d %s','headerLines',1);
%     fclose(fileID);


%     dat = table(C1{1}, C1{2}, [cell2mat(C1(3)) cell2mat(C1(4))], 'VariableNames',{'OriginFileName','ImageName','XY'});
    pid = PIDs_Att_Exp{p};

    p_dat = AttExp(strcmp(AttExp(:,1),pid),:);
%     fnames = tps{p,1};  % odd matrix

    task_dat = cell2mat(p_dat(:,2));
    inames = p_dat(:,3);  % even matrix
    xs = cell2mat(p_dat(:,4)); 
    ys = cell2mat(p_dat(:,5)); 
    
    for t=1:length(tasks)
        task = task_dat == t;
        dat = table(task_dat(task), inames(task), [xs(task) ys(task)], 'VariableNames',{'OriginFileName','ImageName','XY'});

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
            output_name = [out_dir img_names{i} '_Patient_' patients{p} '_observed_touchpoints_' tasks{t} '.mat'];

            save(output_name,'img_dataset');
        end % shape loop
    end % task loop
end %patient loop
