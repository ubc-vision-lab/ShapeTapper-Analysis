clear; clc;

working_dir = pwd; % make sure this directory has /data, /images, and config.txt
demographics_dir = '\demographics';
output_dir = "output_dir";

nn = 1;

%% scrape all subjects in the files
subject_files = getAllFiles(strcat(working_dir, demographics_dir));
subjects = cell(1,numel(subject_files));
num_subjects = 1;
for i = 1:numel(subject_files)
    subject_id = strsplit(subject_files{i},'\');
    subject_id = strsplit(subject_id{end},'_'); % file names always start with <subjectID>_<name>.txt
    subject_id = subject_id{1};
    if ~any(strcmp(subjects,subject_id)) % only add if it's not already there
        subjects{num_subjects} = subject_id;
        num_subjects = num_subjects + 1;
    end
end
subjects(~cellfun('isempty',subjects)); % remove all empty cells

output = {};

for j = 1:numel(subjects) % go through all subjects mined from above
    % potentially we can all every subject we've completed into a mat file
    % so we know what we've already finished.
    
    curr_subject = subjects{j};

    %scrape demo file
    fileID = fopen(strcat(working_dir, demographics_dir, '\', curr_subject, '_demographic.txt'), 'r');
    demoData = textscan(fileID,'%s','Delimiter','\n');
    fclose(fileID);
    demoData = demoData{1}; %header?
    temp = strsplit(demoData{2},','); % values

    screen_w = str2num(temp{5}); %change to screen width of device
    screen_h = str2num(temp{6}); %change to screen height of device
    screen_dpi = str2num(temp{7}); %change to screen DPI of device
    config_file = temp{8}; %get config file for subject

    %% grab all the images
    %grab all images
    fileList = getAllFiles(strcat(working_dir, '\images'));
    img_dict = containers.Map;
    %makes a dictionary containing all the shape files and their dimensions
    for image = 1:numel(fileList)
        image_data=imread(fileList{image});

        [height, width, ~] = size(image_data);
        image_name = strsplit(fileList{image},'\');
        image_name = strsplit(image_name{end},'.');
        image_name = image_name{1};
        img_dict(image_name) = [width, height];
    end

    %% calculate transformation data through the configurations
    %go through the config file, then collect and calc transformation data
    fileID = fopen(strcat(working_dir,'\configs\',config_file));
    allData = textscan(fileID,'%s','Delimiter','\n');
    fclose(fileID);
    allData = allData{1};
    trans_info = {};

    for line = 1:numel(allData)
        specs = strsplit(allData{line}, ',');
        block = str2num(specs{1});
        trial = str2num(specs{2});

        x_pos = str2num(specs{15});
        y_pos = str2num(specs{16});

        if str2num(specs{21})
            safety = str2num(specs{20});
            rotation = str2num(specs{19});
            shape = specs{18};
        elseif str2num(specs{29})
            safety = str2num(specs{28});
            rotation = str2num(specs{34});
            shape = specs{26};
        else %35
            safety = str2num(specs{36});
            rotation = str2num(specs{35});
            shape = specs{34};
        end

        new_screenH = screen_h * safety / 100;
        dim = img_dict(shape);
        diag = sqrt(dim(1)^2 + dim(2)^2);
        toScale = diag / new_screenH; %toScale back up
        diag = diag / toScale;

        %border = screen_dpi * 0.375; %0.375 inches = ~ 1cm (so 5mm border around to screen)
        border = screen_dpi; %0.375 inches = ~ 1cm (so 5mm border around to screen)

        margin = screen_w - diag - border;
        toShift = x_pos / 100 * margin;
        shiftAmountX = toShift - (margin / 2);

        margin = screen_h - diag - border;
        toShift = y_pos / 100 * margin;
        shiftAmountY = toShift - (margin / 2);

        plusX = dim(1)/2;
        plusY = dim(2)/2;

        trans_info{block}{trial} = {toScale, -1 * shiftAmountX, -1 * shiftAmountY, 360-rotation, shape, plusX, plusY};

    end

    %% grab the touches
    %go through all data files and collect the touches
    dataList = getAllFiles(strcat(working_dir, '\data'));
    touchData = {};
    n=1;
    for data = 1:numel(dataList)
        filename = strsplit(dataList{data},'\');
        filename = strsplit(filename{end},'.');
        filename = strsplit(filename{1},'_');
        subject_ID = filename{1};
        block_num = str2num(filename{2});

        if isempty(subject_ID)
        elseif(subject_ID == curr_subject) % this matches the subject we're looking at
            fileID = fopen(dataList{data});
            trialData = textscan(fileID,'%s','Delimiter','\n');
            fclose(fileID);
            trialData = trialData{1};
            for trial = 2:numel(trialData)
                temp = strsplit(trialData{trial},',');
                if str2num(temp{2}) == 0
                    touchData{n} = [block_num, str2double(temp{1}), str2double(temp{4}), str2double(temp{5})]; %block, trial, x, y
                    n = n + 1;
                end
            end
        end
    end

    %% generate touch output
    for touch = 1:numel(touchData)
        block = touchData{touch}(1);
        trial = touchData{touch}(2);
        x = touchData{touch}(3);
        y = touchData{touch}(4);
        toScale = trans_info{block}{trial}{1};
        shiftAmountX = trans_info{block}{trial}{2};
        shiftAmountY = trans_info{block}{trial}{3};
        rotation = trans_info{block}{trial}{4};
        shape_name = trans_info{block}{trial}{5};
        plusX = trans_info{block}{trial}{6};
        plusY = trans_info{block}{trial}{7};

        x = x + shiftAmountX;
        y = y + shiftAmountY;

        x = x - screen_w/2;
        y = y - screen_h/2;
        
        x = x * toScale;
        y = y * toScale;

        x_old = x;

        x = x * cosd(rotation) - y * sind(rotation);
        y = x_old * sind(rotation) + y * cosd(rotation);

        x = x + plusX;
        y = y + plusY;
        
        output{nn} = {subject_id, shape_name, x, y};
        nn = nn + 1;
    end
    
    %% make directory and print outputs
    if ~(exist(output_dir,'file'))
        mkdir output_dir;
    end
    strcat(working_dir,'\',output_dir,'\',curr_subject,'_results.txt')
    fid = fopen(strcat(working_dir,'\',output_dir,'\',curr_subject,'_results.txt'),'w');
    for i = 1:length(output)
        %fprintf(fid, '%s\t%s\t%s\r\n',output{i}{2:4});
        fprintf(fid, [num2str(output{i}{2}) '\t' num2str(output{i}{3}) '\t' num2str(output{i}{4}) '\r\n']);
    %     fprintf(fid, [output{i}{1} '\t' num2str(output{i}{2}) '\t' num2str(output{i}{3}) '\t' num2str(output{i}{4}) '\r\n']);
    end
    fclose(fid);
end

