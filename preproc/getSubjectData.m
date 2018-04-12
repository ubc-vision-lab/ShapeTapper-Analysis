% Given a directory, grabs all files and read them, assuming they are valid
% data output from Unity #useless comment 

function data = getSubjectData(subjectID,directory)
    data = [];
    files = getAllFiles(directory);
    for i = 1:numel(files)
        filename = strsplit(files{i},filesep);
        filename = filename{end};
        if contains(filename,subjectID)
            fileID = fopen(files{i},'r');
            data = [data; textscan(fileID,"%d %d %f %d %d",'HeaderLines',1,'Delimiter',',')];
            fclose(fileID);
        end
    end
end