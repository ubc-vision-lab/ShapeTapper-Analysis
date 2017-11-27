function demographic_data = getDemographicData(filepath)
    fileID = fopen(filepath, 'r');
    header_format_spec = '%s';
    num_header_columns = 8;
    demoHeaders = textscan(fileID, header_format_spec, num_header_columns, 'Delimiter',',');
    demographic_data = textscan(fileID,'%s %s %s %d %d %d %d %s','Delimiter',',');
    fclose(fileID);
end