% read the results file
% create the bins with all the indices
% grab the image
% winning?
clear; clc;

base_dir = pwd; % change this to the directory you want. Everything assumed to be from here.
subj_data_map = containers.Map;
image_files = getAllFiles(fullfile(base_dir,'images'));

results_name = 'results';
output_name = 'output_overlays';
per_trial_images = true;
if per_trial_images
    output_name = strcat(output_name,'_trial');
    results_name = strcat(results_name,'_trial');
end
scaled = true;
absolute = true;
if scaled
    results_name = strcat(results_name,'_scaled');
    output_name = strcat(output_name,'_scaled');
end
if absolute
    results_name = strcat(results_name,'_absolute');
    output_name = strcat(output_name,'_absolute');
end
if(~exist(output_name,'dir'))
    mkdir(output_name);
end
results_dir = fullfile(base_dir,results_name);
output_dir = fullfile(base_dir,output_name);
result_files = getAllFiles(results_dir);
% result_files = {fullfile(results_dir,'test_results_old_rounded.txt')};

%% get touch results

% make folders for the images for each subject
% something like results/subject_images/*

for i = 1:size(result_files)
    result_file = result_files{i};
    if ( exist(result_file,'file') && contains(result_file,'results') ) % it's a result file
        subjectID = split(result_file,filesep);
        subjectID = subjectID(end);
        subjectID = split(subjectID,'_');
        subjectID = subjectID{1}; % this is the actual subject ID
        result_fid = fopen(result_file);
        subj_data_map(subjectID) = textscan(result_fid,'%s %f %f %d %d'); % readtable is slow
        fclose(result_fid);
    end
end

%% Generate the images
for k = keys(subj_data_map)
    subj_touch_data = subj_data_map(k{1});
    %here, we'll assume it's a 1x5 array since we're using textscan
    for i = 1:size(subj_touch_data{1},1) % this gets us the length of the column
        image_name = subj_touch_data{1}(i);
        for j = 1:size(image_files,1)
            if contains(image_files(j),image_name) % the subject ID is in here (assumption, the sID doesn't appear in the directory structure)
                image = imread(image_files{j});
                shown_img = imshow(image_files{j});
                hold on;
                [img_y, img_x, ~] = size(image);
                overlay = plot(subj_touch_data{2}(i), img_y-subj_touch_data{3}(i),'om');
                
                subject_output_dir = fullfile(output_dir,k{1});
                subject_output = fullfile(subject_output_dir,strcat(k{1},'_',num2str(subj_touch_data{4}(i)),'_',num2str(subj_touch_data{5}(i))));
                if(~exist(subject_output_dir,'dir'))
                    mkdir(subject_output_dir);
                end
                saveas(overlay,subject_output,'png');
        %         pcolor(x,y,touch_rounded);
                hold off
                break;
            end
        end
    end
end