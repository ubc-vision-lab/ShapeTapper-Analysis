function results = getSubjectResults(subjectID,directory)
    results = [];
    files = getAllFiles(directory);
    for i = 1:numel(files)
        filename = strsplit(files{i},filesep);
        filename = filename{end};
        if contains(filename,subjectID)
            fileID = fopen(files{i},'r');
            results = [results; textscan(fileID,"%d %d %f %d %d",'HeaderLines',1,'Delimiter',',')];
            fclose(fileID);
        end
    end
end