clear

% patients = {'S01','S02','S03','S04','S06','S07',...
%             'S08','S09','S11','S13','S14','S17',...
%             'S20','S21','S22','S23','S24'};
patients = {'S25','S26','S27',...
            'S28','S29','S30',...
            'S31','S32','S33'};
patient_nums = [25 26 27 28 29 30 31 32 33];
        
in_path  = 'D:\ShapeTapper-Analysis\';
out_path = 'D:\ShapeTapper-Analysis\';

for p=1:length(patients)
    
    fileID = fopen([in_path patients{p} '\' num2str(patient_nums(p)) '_Overall_Results.txt']);
    C1 = textscan(fileID,'%s %f %f %d %s','headerLines',1);
    fclose(fileID);

    img_names = ["blake_01","blake_03","blake_04","blake_06","blake_07","blake_08","blake_09","blake_10","blake_11","blake_12"];

%     dat = table(C1{1}, C1{2}, [cell2mat(C1(3)) cell2mat(C1(4))], 'VariableNames',{'OriginFileName','ImageName','XY'});

    fnames = C1{1};  % odd matrix
    inames = C1{5};  % even matrix
    xs = C1{2}; 
    ys = C1{3}; 
    dat = table(fnames, inames, [xs ys], 'VariableNames',{'OriginFileName','ImageName','XY'});
    
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
        output_name = [out_dir img_names{i} '_Patient_' patients{p} '_observed_touchpoints.mat'];

        save(output_name,'img_dataset');

    end %shape loop
end %patient loop