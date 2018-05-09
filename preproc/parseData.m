patient = 'MC1';

fileID = fopen(['..\' patient '\aggregated_observations\Overall_Results.txt']);
C1 = textscan(fileID,'%s %s %f %f','headerLines',1);
fclose(fileID);

img_names = ["solo3","solo5","solo6","solo7","solo9","solo10","solo11","solo12",...
             "blake_01","blake_03","blake_04","blake_06","blake_07","blake_08","blake_09","blake_10","blake_11","blake_12"];

dat = table(C1{1}, C1{2}, [cell2mat(C1(3)) cell2mat(C1(4))], 'VariableNames',{'OriginFileName','ImageName','XY'});

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
    
    output_name = ['..\' patient '\aggregated_observations\' img_names{i} '_Patient_' patient '_aggregated_observations.mat'];
    
    save(output_name,'img_dataset');
    
end