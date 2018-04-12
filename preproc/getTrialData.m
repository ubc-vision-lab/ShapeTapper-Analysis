function trial_data = getTrialData(dataFile)
    datafID = fopen(dataFile,'r');
    header_format_spec = '%s';
    num_header_columns = 5;
    trial_headers = textscan(datafID, header_format_spec, num_header_columns, 'Delimiter',',');
    trial_data = textscan(fileID,'%d %d %d %d %d','Delimiter',',');
    fclose(fileID);
end