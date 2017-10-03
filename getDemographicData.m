function demographic_data = getDemographicData(filepath)
fileID = fopen(filepath, 'r');
demoData = textscan(fileID,'%s','Delimiter','\n');
fclose(fileID);
demoData = demoData{1}; %header?
demographic_data = strsplit(demoData{2},','); % values
end