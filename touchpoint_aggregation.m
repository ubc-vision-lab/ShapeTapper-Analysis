clear; clc;

subjects = getAllSubjects(strcat(pwd, filesep, 'demographics'));
results_row = 4;
results_col = 5;

for i = 1:numel(subjects)
    demographic_data = getDemographicData(strcat(pwd,filesep, 'demographics',filesep,subjects{i},'_demographic.txt'));
    screen_w = str2num(demographic_data{5}); %change to screen width of device
    screen_h = str2num(demographic_data{6}); %change to screen height of device
    screen_dpi = str2num(demographic_data{7}); %change to screen DPI of device
    config_file = demographic_data{8}; %get config file for subject
    config_name = strsplit(config_file,'.');
    config_name = config_name{1};
    save_dir = strcat(subjects{i},'_',config_name);
    % grab all the results
    results = getSubjectResults(subjects{i}); % each row is a block, and the columns inside correspond to trials
    
    % grab the configs
    configurations = getConfigurations(config_file); % each row is a trial (sorry!), also it's in string (sorry!)
    if ~exist(save_dir,'dir')
        mkdir(save_dir);
    end
    
    % create a bin for every image
    % find the closest shape to the touchpoint (center?)
    image_bins = containers.Map;
    images = getAllFiles(strcat('.',filesep,'images'));
        
    for j = 1:length(configurations) % go through every configuration
        [events, block, trial] = getImageData(configurations{j});
        if block > size(results,1)
            break;
        end
        
        % find the closest event, add to bin
        closest_distance = inf;
        closest_image = '';
        closest_touchpt = [0, 0];
        diag_max = (screen_h * events{i,1}{5}) / 100;
        for k = 1:length(events,1)
            touchpt = [results{block,results_col}(trial),results{block,results_col}(trial)];
            event_center = [events{i}{2}/100 * (screen_w - diag_max) + diag_max/2 , events{i}{3}/100 * (screen_h - diag_max) + diag_max/2];
            distance = norm(touchpt-event_center);
            if(distance < closest_distance)
                closest_distance = distance;
                closest_image = events{i}{1};
                closest_touchpt = touchpt;
            end
        end
        % rotate the touchpoint relative to the rotation of the image
        
        image_bins(closest_image) = [image_bins(closest_image); closest_touchpt];
    end
    
    % go through all the bins and generate the images, with some generous
    % overhead
    [trial_images, trial_alpha] = composeImage(events,strcat(pwd,filesep,'images'),screen_w,screen_h);
        trial_background = imshow(192*ones(screen_h,screen_w,3,'uint8'));
        hold on;
        % add a point for where the touch occurred
        margin = 5;
        
        %add the touchpoint to the picture
        touchpixel = 255 * ones(2*margin+1,2*margin+1,3,'uint8');
        touchpixel(:,:,2) = 0;
        if results{block,results_col}(trial)~=0 || results{block,results_row}(trial) ~= 0
            trial_images((screen_h - results{block,results_col}(trial)-margin):(screen_h - results{block,results_col}(trial)+margin),...
            results{block,results_row}(trial)-margin:results{block,results_row}(trial)+margin,:) = touchpixel; % some colour
        
            %apply also to the alpha or else your touchpoint won't appear!
            trial_alpha((screen_h - results{block,results_col}(trial)-margin):(screen_h - results{block,results_col}(trial)+margin),...
            results{block,results_row}(trial)-margin:results{block,results_row}(trial)+margin) = 255 * ones(2*margin+1,2*margin+1,'uint8');
        end
        
        %show the figure
        trial_view = imshow(trial_images);
        set(trial_view,'AlphaData',trial_alpha);
        save_name = strcat(pwd,filesep,save_dir,filesep,save_dir,'_',num2str(block),'_',num2str(trial));
        saveas(trial_view,save_name,'png')
        hold off;
    % generate the target image
    
    % overlay the touchpoint on that image from the results file
    % wait for user input, or save to disk
end