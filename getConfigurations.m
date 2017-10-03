function config_data = getConfigurations(configName)
    config_data = [];
    files = getAllFiles(strcat(pwd, filesep,'configs'));
    for i = 1:numel(files)
        filename = strsplit(files{i},filesep);
        filename = filename{end};
        if contains(filename,configName)
            fileID = fopen(files{i},'r');
            config_data = textscan(fileID,"%s",'Delimiter','\n');
            fclose(fileID);
            config_data = config_data{1};
            for j = 1:numel(config_data)
                config_data{j} = strsplit(config_data{j},',');
            end
        end
    end
end