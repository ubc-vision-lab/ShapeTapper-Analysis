function p_vals = typeI_est_rsdt(x_beta_pdf, y_beta_pdf, n, n_trials)
%typeI_est_zscores Estimate Type I Error rates for t-diff for given beta PDFs

p_vals = zeros(1,n_trials);

for i=1:n_trials
    sample_x = random(x_beta_pdf, [1, n+1]);
    sample_y = random(y_beta_pdf, [1, n+1]);
    [~, p_vals(i)] = rsdt_crawford(sample_x(1), sample_x(2:end), sample_y(1), sample_y(2:end));
end
    
end