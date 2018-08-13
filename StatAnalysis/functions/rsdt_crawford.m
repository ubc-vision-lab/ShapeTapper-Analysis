function [t_diff, p_val, r] = rsdt_crawford(patient_task_x, sample_task_x, patient_task_y, sample_task_y)
%tdiff_crawford Get t-score diff + pval according to Crawford & Howell 1998

n = length(sample_task_x);

corr = corrcoef(sample_task_x, sample_task_y,'Rows','complete');
r = corr(1,2);

Z_x = get_z_score(patient_task_x, sample_task_x);
Z_y = get_z_score(patient_task_y, sample_task_y);

a = (1-r)*(1-r^2);
b = (1-r) * ( (4*(n-1)^2) + 4*(1+r)*(n-1) + (1+r)*(5+r) );

c = -2 * (Z_x - Z_y)^2 * ((n*(n-1)^2)/(n+1));

t_diff = -sqrt( (-b + sqrt(b^2 - 4*a*c)) / (2*a) );
p_val = tcdf(t_diff,n);
            
end

