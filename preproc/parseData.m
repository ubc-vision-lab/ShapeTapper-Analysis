fileID = fopen("..\MC\aggregated_observations\Overall_Results.txt");
C = textscan(fileID,'%s %s %f %f','headerLines',1);
fclose(fileID);

img_names = ["solo3","solo5","solo6","solo7","solo9","solo10","solo11","solo12",...
             "blake_01","blake_03","blake_04","blake_06","blake_07","blake_08","blake_09","blake_10","blake_11","blake_12"];

dat = table(C{1}, C{2}, [cell2mat(C(3)) cell2mat(C(4))], 'VariableNames',{'OriginFileName','ImageName','XY'});

%dat = dat(~strcmp(dat.OriginFileName, 's_solo_blake_10b36p_3m_Y8hZ'),:);

for i=1:length(img_names)
    
    shape_dat = dat(strcmp(dat.ImageName, img_names(i)),:);
    img_dataset = shape_dat.XY;
    
    output_name = ['..\MC\aggregated_observations\' img_names{i} '_Patient_MC_aggregated_observations.mat'];
    
    save(output_name,'img_dataset');
    
end