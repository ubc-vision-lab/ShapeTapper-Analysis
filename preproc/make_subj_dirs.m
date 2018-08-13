clear

patients = {'S06',...
            'S07','S08','S09','S10','S11','S12',...
            'S13','S14','S15','S16','S17','S18',...
            'S19','S20','S21','S22','S23','S24'};

% in_path  = 'D:\ShapeTapper-Analysis\';
out_path = 'E:\ShapeTapper-Analysis\';

for p=1:length(patients)
    
    out_dir = [out_path patients{p} '\observed_touchpoints\'];

    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    out_dir = [out_path patients{p} '\uniform_points\'];

    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    out_dir = [out_path patients{p} '\spatial_analysis\'];

    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

end %patient loop