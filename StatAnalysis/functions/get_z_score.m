function z_score = get_z_score(x,sample)
%get_z_score Get z-score based on mean and sample stdev

mean = nanmean(sample);
stdev = nanstd(sample,0);
z_score = (x-mean)/stdev;
end

