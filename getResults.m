function results = getResults(subjectID,results_dir,mode)
    success = false;
    results_files = getAllFiles(results_dir);
    for i = 1:numel(results_files)
        if ( contains(results_files{i},subjectID) && contains(results_files{i},'results') ) % it's a result file
            if mode % table mode
                results = readtable(results_files{i},'ReadVariableNames',false);
                results.Properties.VariableNames = {'image','x','y','block','trial'};
                success = true;
            else % cell mode
                result_fid = fopen(results_files{i});
                results = textscan(result_fid,'%s %f %f %d %d'); % readtable is slow
                success = true;
                fclose(result_fid);
            end
        end
    end
    if ~success
        results = {};
    end
end