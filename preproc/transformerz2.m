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
absolute = true;
output_dir = 'results';
if per_trial_transformation
    output_dir = strcat(output_dir,'_trial');
end
if scaled
    output_dir = strcat(output_dir,'_scaled');
end
if absolute
    output_dir = strcat(output_dir,'_absolute');
end
output_dir = fullfile(pwd,output_dir);

nn = 1;

%% scrape all subjects in the files
subjects = getAllSubjects(demographics_dir);

% % grab subject IDs from a file
% subjects = textscan(fopen(fullfile(pwd,'Patient_2_IDs.txt')),'%s');
% subjects = subjects{1};

%% go through all subjects and grab the demographic data
for j = 1:numel(subjects) % go through all subjects mined from above
    % potentially we can all every subject we've completed into a mat file
    % so we know what we've already finished.
    output = {};
    nn = 1;
    
    curr_subject = subjects{j};
    if strcmp(curr_subject,'test')
        nn = nn;
    end

    %scrape demo file
    fileID = fopen(fullfile(demographics_dir, strcat(curr_subject, '_demographic.txt')), 'r');
    demoData = textscan(fileID,'%s','Delimiter','\n');
    fclose(fileID);
    demoData = demoData{1}; % this gets all the rows
    temp = strsplit(demoData{2},','); % values

    screen_w = str2double(temp{5}); %change to screen width of device
    screen_h = str2double(temp{6}); %change to screen height of device
    screen_dpi = str2double(temp{7}); %change to screen DPI of device
    config_file = temp{8}; %get config file for subject
    border = 0;

    %% get the touches
    %go through all data files and collect the touches
    dataList = getAllFiles(data_dir);
    touchData = {};
    touchData_map = containers.Map;
    n=1;
    for data = 1:numel(dataList) % All data files in the folder
        filename = strsplit(dataList{data},filesep); % data files are formatted as '<subject_ID>_<block_num>.txt' (no angle braces)
        filename = strsplit(filename{end},'.');
        filename = strsplit(filename{1},'_');
        subject_ID = filename{1};
        block_num = str2double(filename{2});

        if isempty(subject_ID) % the filename started with an underscore?
            continue
        elseif(subject_ID == curr_subject) % this matches the subject we're looking at
            fileID = fopen(dataList{data});
            % this section can likely be replaced by readtable()
            trialData = textscan(fileID,'%s','Delimiter','\n'); % erm. by line? dunno why
            fclose(fileID);
            trialData = trialData{1}; % it's a layer deeper because of textscan so that we can split it
            for trial = 2:numel(trialData) % actual data rows, first row is headers
                temp = strsplit(trialData{trial},','); % split the data
                if str2double(temp{2}) == 0
                    % touch points from Unity are from bottom left corner
                    % (Input.GetTouch(0).position)
                    touchData{n} = [block_num, str2double(temp{1}), str2double(temp{4}), str2double(temp{5})]; %block, trial, x, y
                    touchData_map(strcat(filename{2}, '+', temp{1})) = [str2double(temp{4}), str2double(temp{5})];
                    n = n + 1;
                end
            end
        end
    end
    
    %% grab all the images
    %grab all images
    fileList = getAllFiles(images_dir);
    img_dict = containers.Map;
    %makes a dictionary containing all the shape files and their dimensions
    for image = 1:numel(fileList)
        image_data=imread(fileList{image});

        [height, width, ~] = size(image_data);
        image_name = strsplit(fileList{image},filesep);
        image_name = strsplit(image_name{end},'.');
        image_name = image_name{1};
        img_dict(image_name) = [width, height];
    end

    %% calculate transformation data through the configurations
    %go through the config file, then collect and calc transformation data
    fileID = fopen(strcat(config_dir,filesep,config_file));
    allData = textscan(fileID,'%s','Delimiter','\n');
    fclose(fileID);
    allData = allData{1};
    trans_info = {};

    % calculate the info for the images
    % trans_info is the touch transformation data
    for line = 1:numel(allData)
        specs = strsplit(allData{line}, ',');
        block = str2double(specs{1});
        trial = str2double(specs{2});
        
        event_info = {};
        if(~isempty(specs{18}))
            shape = specs{18};
            rotation = str2double(specs{19});
            safety = str2double(specs{20});
            img_dim = img_dict(shape);
            event_info{1} = calculate_event_info(screen_w, screen_h, screen_dpi, safety, shape, img_dim, border, str2double(specs{15}), str2double(specs{16}), rotation );
        end
        if(~isempty(specs{26}))
            shape = specs{26};
            rotation = str2double(specs{27});
            safety = str2double(specs{28});
            img_dim = img_dict(shape);
            event_info{2} = calculate_event_info(screen_w, screen_h, screen_dpi, safety, shape, img_dim, border, str2double(specs{42}), str2double(specs{43}), rotation );
        end
        if(~isempty(specs{34}))
            shape = specs{34};
            rotation = str2double(specs{35});
            safety = str2double(specs{36});
            img_dim = img_dict(shape);
            event_info{3} = calculate_event_info(screen_w, screen_h, screen_dpi, safety, shape, img_dim, border, str2double(specs{44}), str2double(specs{45}), rotation );
        end
        trans_info{block}{trial} = event_info;

    end
   

    %% generate touch output
    for touch = 1:numel(touchData)
        block = touchData{touch}(1);
        trial = touchData{touch}(2);
        
        touchpoint = [touchData{touch}(3), touchData{touch}(4)];
        image_relative_touchpoint = [0 0];
        closest_distance = inf;
        shape_name = '';
        
%         touch = touchData_map(strcat(specs{1},'+',specs{2}));
%         for event_position_index = 1:size(event_info,1)
%         end
        
        % TODO: Change this to select the closest shape -- it's a deeper
        % cell array than before!
        for event = 1:size(trans_info{block}{trial},2)
            toScale = trans_info{block}{trial}{event}{1};
            shift = [trans_info{block}{trial}{event}{2} trans_info{block}{trial}{event}{3}];
            rotation = trans_info{block}{trial}{event}{4};
            plusX = trans_info{block}{trial}{event}{6};
            plusY = trans_info{block}{trial}{event}{7};
            
            difference_vector = touchpoint + shift; %put touchpoint into image relative position
%             difference_vector = difference_vector - [screen_w/2, screen_h/2]; % old math
            difference_vector = difference_vector * toScale; % scale it to the spot relative to the original image
            distance = norm(difference_vector);
            if(distance < sqrt(plusX^2 + plusY^2) && distance < closest_distance)
                closest_distance = distance;
                shape_name = trans_info{block}{trial}{event}{5};
                R = [cosd(rotation) sind(rotation); -sind(rotation) cosd(rotation)];
                image_relative_touchpoint = difference_vector * R;
                image_relative_touchpoint = image_relative_touchpoint + [plusX, plusY];
            end
        end
        if ~isempty(shape_name)
            if(per_trial_transformation)
                output{nn} = {curr_subject, shape_name, image_relative_touchpoint(1), image_relative_touchpoint(2), block, trial};
            else
                output{nn} = {curr_subject, shape_name, image_relative_touchpoint(1), image_relative_touchpoint(2)};
            end
            nn = nn + 1;
        end
    end
    
    %% make directory and print outputs
    if ~(exist(output_dir,'file'))
        mkdir(output_dir);
    end
    fid = fopen(strcat(output_dir,filesep,curr_subject,'_results.txt'),'w');
    for i = 1:length(output)
        %fprintf(fid, '%s\t%s\t%s\r\n',output{i}{2:4});
        if(per_trial_transformation)
            fprintf(fid, [num2str(output{i}{2}) '\t' num2str(output{i}{3}) '\t' num2str(output{i}{4}) '\t' num2str(output{i}{5}) '\t' num2str(output{i}{6}) '\r\n']);
        else
            fprintf(fid, [num2str(output{i}{2}) '\t' num2str(output{i}{3}) '\t' num2str(output{i}{4}) '\r\n']);
        end
    %     fprintf(fid, [output{i}{1} '\t' num2str(output{i}{2}) '\t' num2str(output{i}{3}) '\t' num2str(output{i}{4}) '\r\n']);
    end
    fclose(fid);
end
