function demographic_data = getDemographicData(filepath)
    fileID = fopen(filepath, 'r');
    demoData = textscan(fileID,'%s','Delimiter','\n');
    fclose(fileID);
    demoData = demoData{1}; % all the rows are in cell, get that one cell
    demographic_data = strsplit(demoData{2},','); % values
end