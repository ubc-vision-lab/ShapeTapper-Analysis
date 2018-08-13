function [t_diff, p_val] = tdiff_crawford(patient_task_x, sample_task_x, patient_task_y, sample_task_y)
%tdiff_crawford Get t-score diff + pval according to Crawford & Howell 1998

n = length(sample_task_x);

Z_x = get_z_score(patient_task_x, sample_task_x);
Z_y = get_z_score(patient_task_y, sample_task_y);
r_xy = corrcoef(sample_task_x, sample_task_y,'Rows','complete');

t_diff = (Z_x - Z_y) / sqrt( (2-2*r_xy(1,2)) * ((n+1)/n));
p_val = tcdf(t_diff,n);
            
end

