clear; clc;

% pwd is the Present Working Directory (the folder you're in)
% filesep is the OS agnostic separator
% (UNIX-like (OSX, Linux) '/', Windows '\')
working_directory = pwd; % change this to whichever directory you want to work on
demographics_dir = fullfile(working_directory,'demographics');
images_dir = fullfile(working_directory,'images');
data_dir = fullfile(working_directory,'data_new');
config_dir = fullfile(working_directory,'configs');
per_trial_transformation = true;
scaled = true;
output_dir = 'results';
if per_trial_transformation
    output_dir = strcat(output_dir,'_trial');
end
if scaled
    output_dir = strcat(output_dir,'_scaled');
end
output_dir = fullfile(pwd,output_dir);

%% grab all the demographic data
% grab the subject data, allow user to change set of subjects here.

demographics_data = cellfun(@getDemographicData,getAllFiles(demographics_dir));
% subjectIDs = cellfun(@(x) x{1}(1),demographics_data);
% Uncomment next line if you only want certain subjects processed
% subjectIDs = {'4rD9', 'hQx3'} % put the subject IDs you want here
subjectIDs = textscan(fopen(fullfile(pwd,'Patient_2_IDs.txt')),'%s');
subjectIDs = subjectIDs{1};

%% open all data files, grab all data for this subject
touchPoint_image_map = container.Map;

data_filenames = getAllFiles(data_dir);
for i = 1:numel(data_filenames)
    
    if(any(contains(data_files,subjectIDs)))
        
    end
end

