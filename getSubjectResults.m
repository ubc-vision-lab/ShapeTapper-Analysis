function results = getSubjectResults(subjectID)
    results = [];
    files = getAllFiles(strcat(pwd, '\data'));
    for i = 1:numel(files)
        filename = strsplit(files{i},'\');
        filename = filename{end};
        if contains(filename,subjectID)
            fileID = fopen(files{i},'r');
            results = [results; textscan(fileID,"%u16 %u16 %f %u16 %u16",'HeaderLines',1,'Delimiter',',')];
            fclose(fileID);
        end
    end
end