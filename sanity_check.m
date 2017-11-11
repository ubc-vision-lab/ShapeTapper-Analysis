clear; clc;

subjects = getAllSubjects(strcat(pwd, filesep, 'demographics'));
subjects = {'Dgwm'};
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
    for j = 1:length(configurations)
        [events, block, trial] = getImageData(configurations{j});
        if block > size(results,1)
            break;
        end
        % we've already made this image, next image
%         if exist(strcat(pwd,filesep,save_dir,filesep,save_dir,'_',num2str(block),'_',num2str(trial),'.png'),'file')
%             continue;
%         end
        [trial_images, trial_alpha] = composeImage(events,strcat(pwd,filesep,'images'),screen_w,screen_h, screen_dpi);
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
    end
    % generate the target image
    
    % overlay the touchpoint on that image from the results file
    % wait for user input, or save to disk
end