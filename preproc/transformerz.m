working_dir = 'C:\Users\Arnold\Desktop\UBC Vision Lab\d usb dump\rob stuff\shapetapper data analysis\iOS DATA'; % make sure this directory has /data, /images, and config.txt

nn = 1;
output = {};
for user = 1:34 %1:length(user_configs) %it's a mat file. Don't know why.

    subject_id = user_configs{user,1};
    config_file = user_configs{user,2}; % corresponding configuration


    %scrape demo file
    fileID = fopen(strcat(working_dir, '\data', '\', subject_id, '_demographic.txt'));
    demoData = textscan(fileID,'%s','Delimiter','\n');
    fclose(fileID);
    demoData = demoData{1}; %header?
    temp = strsplit(demoData{2},',');

    screen_w = str2num(temp{5}); %change to screen width of device
    screen_h = str2num(temp{6}); %change to screen height of device
    screen_dpi = str2num(temp{7}); %change to screen DPI of device

    fileList = getAllFiles(strcat(working_dir, '\images'));
    dict = containers.Map;

    %makes a dictionary containing all the shape files and their dimensions
    for images = 1:numel(fileList)
        image2=imread(fileList{images});

        [height, width, ~] = size(image2);
        C = strsplit(fileList{images},'\');
        C = strsplit(C{end},'.');
        C = C{1};
        dict(C) = [width, height];
    end


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

        x_pos = str2num(specs{11});
        y_pos = str2num(specs{12});

        if str2num(specs{17})
            safety = str2num(specs{16});
            rotation = str2num(specs{15});
            shape = specs{14};
        elseif str2num(specs{26})
            safety = str2num(specs{25});
            rotation = str2num(specs{24});
            shape = specs{23};
        else %35
            safety = str2num(specs{34});
            rotation = str2num(specs{33});
            shape = specs{32};
        end

        new_screenH = screen_h * safety / 100;
        dim = dict(shape);
        diag = sqrt(dim(1)^2 + dim(2)^2);
        toScale = diag / new_screenH; %toScale back up
        diag = diag / toScale;

        border = screen_dpi * 0.375; %0.375 inches = ~ 1cm (so 5mm border around to screen)

        margin = screen_w - diag - border;
        toShift = x_pos / 100 * margin;
        shiftAmountX = toShift - (margin / 2);

        margin = screen_h - diag - border;
        toShift = y_pos / 100 * margin;
        shiftAmountY = toShift - (margin / 2);

        plusX = dim(1)/2;
        plusY = dim(2)/2;

        trans_info{block}{trial} = {toScale, -1 * shiftAmountX, -1 * shiftAmountY, 360-rotation, shape, plusX, plusY}

    end

    %go through all data files and collect the touches
    dataList = getAllFiles(strcat(working_dir, '\data'));
    touchData = {};
    n=1;
    for data = 1:numel(dataList)
        C = strsplit(dataList{data},'\');
        C = strsplit(C{end},'.');
        C = strsplit(C{1},'_');
        sub_id = C{1};
        C = str2num(C{2});

        if isempty(C)

        elseif(sub_id == subject_id)  
            fileID = fopen(dataList{data});
            trialData = textscan(fileID,'%s','Delimiter','\n');
            fclose(fileID);
            trialData = trialData{1};
            for trial = 2:numel(trialData)
                temp = strsplit(trialData{trial},',');
                if str2num(temp{2}) == 0
                    touchData{n} = [C, str2num(temp{1}), str2num(temp{4}), str2num(temp{5})]; %block, trial, x, y
                    n = n + 1;
                end    
            end
        end
    end

    %generate touch output
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
end


fid = fopen(['all_results.txt'],'w');
for i = 1:length(output)
    fprintf(fid, [output{i}{1} '\t' num2str(output{i}{2}) '\t' num2str(output{i}{3}) '\t' num2str(output{i}{4}) '\r\n']);
end
fclose(fid);