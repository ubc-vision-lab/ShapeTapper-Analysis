function subjectID = getSubjectID(filename)
    subjectID = strsplit(filename,filesep);
    subjectID = strsplit(subjectID{end},'_'); % file names always start with <subjectID>_<name>[_<etc>].txt
    subjectID = subjectID{1};
end