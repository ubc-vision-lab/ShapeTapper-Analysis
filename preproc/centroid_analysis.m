% bin the results according to image
% 
% grab the centroid for the image
% create the random set using square
% create a random set with normal distr, mu{x,y} = touchpoint_avg{x,y}
% for that image variance{x,y} = touchpoints_var{x,y}
clear; clc

% assumption: All subjectIDs being parsed are from one person
working_dir = pwd; % you can see this to whatever you want -- pwd is the directory you're running the script from)
data_dir = fullfile(working_dir,'data_new');
img_dir = fullfile(working_dir,'images');
output_mat_data = fullfile(working_dir,'Patient_DF_generated_results');
image_space_touchpt_dir = fullfile(working_dir,'results_trial_scaled_absolute');
num_generated_datasets = 10000;
% subjects = getAllSubjects(data_dir);
% subjects = {'5fvs','Dgwm','MWmz','Sbr7','vPzW','Y9uX'}; % for custom list
subjects = textscan(fopen(fullfile(pwd,'Patient_1_IDs.txt')),'%s'); % read from file
subjects = subjects{1};

% get the image space results
if(numel(subjects)>1)
    results = cellfun(@(x) getResults(x,image_space_touchpt_dir,1),subjects,'UniformOutput',false);
else
    results = getResults(subjects{1},image_space_touchpt_dir,1);
end

% If there's only one subject, there's only going to be table, so no need
% to 
if size(results,1) > 1 && size(results,2) == 1
    patient_touch_aggregate = results{1};
    for i = 2:size(results,1)
       patient_touch_aggregate = union(patient_touch_aggregate,results{i});
    end
else % so, in the case the output is a table, the flow ends up here anyway
    patient_touch_aggregate = results;
end

[groups, images] = findgroups(patient_touch_aggregate.image); % grab the groups and corresponding images
mystats = @(x)[numel(x)/2 mean(x) var(x) std(x)];
xy_stats = splitapply(mystats,[patient_touch_aggregate.x,patient_touch_aggregate.y],groups); % xy stats for all images
% columns for xy_stats {'num_touches','mean_x','mean_y','var_x','var_y','std_x','std_y'};

% this lets you extract all the touchpoints for an image
% blake_points = patient_touch_aggregate(strcmp(patient_touch_aggregate.image,'blake_01'),:);

%% grab all images and get their dimensions (and radius
all_img_datasets = containers.Map;

for i = 1:numel(images)
%     generate sets of random data
    tmp_filename = fullfile(img_dir,strcat(images{i},'.png'));
    if exist(tmp_filename,'file')
        tmp_img = imread(tmp_filename);
        img_diag = sqrt(size(tmp_img,2)^2+size(tmp_img,1)^2);
        img_dim_half = [size(tmp_img,2), size(tmp_img,1)]/2; % y is first in imread
    else
        continue; % next loop iter
    end
    
    size_per_dataset = xy_stats(i,1); % number of touchpoints from subject for image
    dataset_mu = xy_stats(i,2:3);
    dataset_var = xy_stats(i,4:5);
    dataset_std = xy_stats(i,6:7);
    img_datasets = zeros(num_generated_datasets,size_per_dataset,2);
    parfor j = 1:num_generated_datasets
        img_dataset = zeros(size_per_dataset,2);
        for k = 1:size_per_dataset % number of touchpoints for subject
            
            % generate a random touchpoint
            new_touchpoint = mvnrnd(dataset_mu,dataset_var);
            while(norm(new_touchpoint - img_dim_half) > img_diag)
                new_touchpoint = rand([dataset_mu-3*,dataset_var));
            end
            img_dataset(k,:) = new_touchpoint;
            % check if it's in our range
            % add it to the point and increment count, otherwise continue
        end
        img_datasets(j,:,:) = img_dataset;
    end
    save(strcat(images{i},'_Patient_DF_generated_outputs'),'img_datasets');
%     all_img_datasets(images{i}) = img_datasets;
end

% save(output_mat_data,'all_img_datasets');



% 100k random seeds for each image
% aggregate all the touchpoints