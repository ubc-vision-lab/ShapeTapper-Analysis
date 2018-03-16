%%%%%%%%%%
% getAllUIDs
% 
% Given a data directory, first part of the filename before the underscore
% '_' character.
% In our case, that should be the unique ID.
% 
%%%%%%%%%%

function all_subs = getAllUIDs(directory)
    subject_files = getAllFiles(directory);
    all_subjects = cell(1,numel(subject_files));
    num_subjects = 1;
    for i = 1:numel(subject_files)
        subject_id = getSubjectID(subject_files{i});
        if ~any(strcmp(all_subjects,subject_id)) % only add if it's not already there
            all_subjects{num_subjects} = subject_id;
            num_subjects = num_subjects + 1;
        end
    end
    all_subs = all_subjects(~cellfun('isempty',all_subjects)); % remove all empty cells
end