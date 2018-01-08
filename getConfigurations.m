function [config_data, config_table] = getConfigurations(configName)
    config_data = [];
    config_table = [];
    files = getAllFiles(strcat(pwd, filesep,'configs'));
    for i = 1:numel(files)
        filename = strsplit(files{i},filesep);
        filename = filename{end};
        if contains(filename,configName)
%             config_table = readtable(filename);
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