clear; clc;

subjectID = "hQx3";

data_dir = pwd;
config_dir = pwd;

demographic_file = subjectID + "_demographic.txt";
demo_file_ID = fopen(data_dir + "\" + demographic_file);
demographic_num_cols = 8;
subject_keys = textscan(demo_file_ID,"%s",demographic_num_cols,'Delimiter',',');
subject_info = textscan(demo_file_ID,"%s %s %s %u %u %u %u %s",'Delimiter',',');
fclose(demo_file_ID);
subject_keys{1}
subject_map = containers.Map(subject_keys{1},subject_info)
subject_map.keys
subject_map.values

%config_data = subject_info(Demographic.Config);
%config_data_ID = fopen(config_data);