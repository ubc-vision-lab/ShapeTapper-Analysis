% 17-11-22 Creates a set of images for each subjectID, where each image is
% an event image overlaid with all touchpoints for that image.

% read the results file
% create the bins with all the indices
% grab the image
% winning?
clear; clc;

results_dir = fullfile(pwd,'results_trial_scaled_absolute');
demographics_dir = fullfile(pwd,'demographics');
result_files = getAllFiles(results_dir);
subj_data_map = containers.Map;
output_name = 'output_overlays_rescaled';
if(~exist(output_name,'dir'))
    mkdir(output_name);
end
output_dir = fullfile(pwd,output_name);

image_files = getAllFiles(fullfile(pwd,'images'));

%% get touch results

% make folders for the images for each subject
% something like results/subject_images/*

for i = 1:size(result_files)
    result_file = result_files{i};
    if ( exist(result_file,'file') && contains(result_file,'results.txt') ) % it's a result file
        subjectID = split(result_file,filesep);
        subjectID = subjectID(end);
        subjectID = split(subjectID,'_');
        subjectID = subjectID{1}; % this is the actual subject ID
        result_fid = fopen(result_file);
        result = textscan(result_fid,'%s %f %f %d %d'); % readtable is slow
        subj_data_map(subjectID) = result;
        fclose(result_fid);
    end
end

%% sort the results by image

for k = keys(subj_data_map)
    subj_touch_data = subj_data_map(k{1});
    subj_touch_data_srtd = containers.Map;
    for i = 1:size(subj_touch_data{1},1)
        if(subj_touch_data_srtd.isKey(subj_touch_data{1}{i}))
            subj_touch_data_srtd(subj_touch_data{1}{i}) = [subj_touch_data_srtd(subj_touch_data{1}{i}); subj_touch_data{2}(i) subj_touch_data{3}(i)];
        else
            subj_touch_data_srtd(subj_touch_data{1}{i}) = [subj_touch_data{2}(i) subj_touch_data{3}(i)];
        end
    end
    
    % now that each set of images has its own touchpoints, start making
    % images out of them, and use a heatmap to map them (round!)
    
    for imgkey = keys(subj_touch_data_srtd)
        img_x = 0;
        img_y = 0;
        for j = 1:size(image_files,1)
            if strfind(image_files{j},imgkey{1})
                image = imread(image_files{j});
                shown_img = imshow(image_files{j});
                hold on;
                 [img_y, img_x, ~] = size(image);
            end
        end
        
        
        touch_rounded = round(subj_touch_data_srtd(imgkey{1}));
        overlay = plot(touch_rounded(:,1), img_y-touch_rounded(:,2),'+m');
        
%         n = min(img_x,img_y);
%         xi = linspace(0,img_x,n);
%         yi = linspace(0,img_y,n);
%         xr = interp1(xi,1:numel(xi),touch_rounded(:,1),'nearest')';
%         yr = interp1(yi,1:numel(yi),touch_rounded(:,2),'nearest')';
%         z = accumarray([xr yr], 1, [n n]);
%         
%         figure(2)
%         surf(z)
        
        subject_output_dir = fullfile(output_dir,k{1});
        subject_output = fullfile(subject_output_dir,strcat(k{1},'_',imgkey{1}));
        if(~exist(subject_output_dir,'dir'))
            mkdir(subject_output_dir);
        end
        saveas(overlay,subject_output,'png');
%         pcolor(x,y,touch_rounded);
        hold off
    end
end