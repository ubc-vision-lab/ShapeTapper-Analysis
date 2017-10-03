clear; clc;

subjects = getAllSubjects(strcat(pwd, filesep, 'demographics'))

for i = 1:1 %numel(subjects)
    demographic_data = getDemographicData(strcat(pwd,filesep, 'demographics',filesep,subjects{i},'_demographic.txt'));
    screen_w = str2num(demographic_data{5}); %change to screen width of device
    screen_h = str2num(demographic_data{6}); %change to screen height of device
    screen_dpi = str2num(demographic_data{7}); %change to screen DPI of device
    config_file = demographic_data{8}; %get config file for subject
    % grab all the results
    results = getSubjectResults(subjects{i}); % each row is a block, and the columns inside correspond to trials
    
    % grab the configs
    configurations = getConfigurations(config_file); % each row is a trial (sorry!), also it's in string (sorry!)
    for j = 1:length(configurations)
        [events, block, trial] = getImageData(configurations{i});
        composeImage(events,strcat(pwd,filesep,'images'),screen_w,screen_h);
    end
    % generate the target image
    
    % overlay the touchpoint on that image from the results file
    % wait for user input, or save to disk
end