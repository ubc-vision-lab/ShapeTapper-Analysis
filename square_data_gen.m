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
output_mat_aggregate_base = '_Patient_DF_aggregated_observations';
output_mat_name_base = '_Patient_DF_generated_results';
output_mat_dir = 'Patient_DF';
if ~exist(output_mat_dir,'dir')
    mkdir(output_mat_dir)
end
image_space_touchpt_dir = fullfile(working_dir,'results_trial_scaled_absolute');
num_generated_datasets = 10000;
rng default % to make this non-repeatable use 'shuffle'
% subjects = getAllSubjects(data_dir);
% subjects = {'5fvs','Dgwm','MWmz','Sbr7','vPzW','Y9uX'}; % for custom list
subjects = textscan(fopen(fullfile(pwd,'Patient_DF_IDs.txt')),'%s'); % read from file
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
%% save the aggregate files into directories
for i = 1:numel(images)
    % save each group
    img_saved_table = patient_touch_aggregate(strcmp(patient_touch_aggregate.image,cellstr(images{i})),2:3);
    img_dataset = table2array(img_saved_table);
    save(fullfile(working_dir,output_mat_dir,strcat(images{i},output_mat_aggregate_base)),'img_dataset');
end

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
    
    number_points_circle = 17;
    angles = linspace(1,360,number_points_circle);
    circlepoly = zeros(number_points_circle,2);
    circlepoly_x = arrayfun(@(x) img_diag*cosd(x)/2+img_dim_half(1),angles);
    circlepoly_y = arrayfun(@(x) img_diag*sind(x)/2+img_dim_half(2),angles);
    
    size_per_dataset = xy_stats(i,1); % number of touchpoints from subject for image
    dataset_mu_x = xy_stats(i,2);
    dataset_mu_y = xy_stats(i,3);
    dataset_var = xy_stats(i,4:5);
    dataset_std_x = xy_stats(i,6);
    dataset_std_y = xy_stats(i,7);
    img_datasets = zeros(num_generated_datasets,size_per_dataset,2);
    sigma = max(ceil(size(tmp_img,2)/dataset_std_x),ceil(size(tmp_img,1)/dataset_std_y))/2; % size of the square distribution
 
    parfor j = 1:num_generated_datasets
        
        img_dataset_x = randi([round(dataset_mu_x-sigma*dataset_std_x) round(dataset_mu_x+sigma*dataset_std_x)],4*sigma*size_per_dataset,1);
        img_dataset_y = randi([round(dataset_mu_y-sigma*dataset_std_y) round(dataset_mu_y+sigma*dataset_std_y)],4*sigma*size_per_dataset,1);
        
        in = inpolygon(img_dataset_x,img_dataset_y,circlepoly_x,circlepoly_y);
        img_dataset = [img_dataset_x(in) img_dataset_y(in)];
        
        if size(img_dataset,1) < size_per_dataset
            "AHHHH"
        end
        img_datasets(j,:,:) = img_dataset(1:size_per_dataset,:);
    end
    save(fullfile(working_dir,output_mat_dir,strcat(images{i},output_mat_name_base)),'img_datasets');
%     all_img_datasets(images{i}) = img_datasets;
end

% save(output_mat_data,'all_img_datasets');



% 100k random seeds for each image
% aggregate all the touchpoints