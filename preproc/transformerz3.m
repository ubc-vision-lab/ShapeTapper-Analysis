% pwd is the Present Working Directory (the folder you're in)
% filesep is the OS agnostic separator
% (UNIX-like (OSX, Linux) '/', Windows '\')
if(~exist('demographics_dir', 'var'))
    demographics_dir = uigetdir;
end
if(~exist('demographics_dir', 'var'))
    images_dir = uigetdir;
end
if(~exist('demographics_dir', 'var'))
    data_dir = uigetdir;
end
if(~exist('demographics_dir', 'var'))
    config_dir = uigetdir;
end
if(~exist('demographics_dir', 'var'))
    output_dir = uigetdir;
end

per_trial_transformation = true;
scaled = true;

%% grab all the demographic data
% grab the subject data, allow user to change set of subjects here.

demographics_data = cellfun(@getDemographicData,getAllFiles(demographics_dir),'UniformOutput',false);
% This next line grabs all subjectIDs
% subjectIDs = cellfun(@(x) x{1}(1),demographics_data);
% Uncomment next line if you only want certain subjects processed
% subjectIDs = {'4rD9', 'hQx3'} % put the subject IDs you want here
% if you have a list of subjects 
subjectIDs = textscan(fopen(fullfile(pwd,'Patient_DF_IDs.txt')),'%s');
subjectIDs = subjectIDs{1};

%% open all data files, grab all data for this subject
touchPoint_image_map = containers.Map;

data_filenames = getAllFiles(data_dir);
for i = 1:numel(data_filenames)
    
    if(any(contains(data_files,subjectIDs)))
        
    end
end

