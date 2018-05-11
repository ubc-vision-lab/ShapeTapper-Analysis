function [obs_dplus, obs_dplus_r, obs_dminus, obs_dminus_r, gen_dplus, gen_dminus] = parse_cdf_struct(cdf_struct)
%PARSECDFSTRUCT Summary of this function goes here
%   Detailed explanation goes here
    unif_cdf = cdf_struct{1,1}{1,1};

    obs_cdf  = cdf_struct{1,1}{1,2};
    obs_cdf_obs_diff = obs_cdf - unif_cdf;
    [obs_dplus, obs_dplus_r]   = max(obs_cdf_obs_diff);
    [obs_dminus, obs_dminus_r] = min(obs_cdf_obs_diff);
    
    gen_cdf = cdf_struct{1,1}{1,3};
    gen_cdf_diff = gen_cdf - unif_cdf;
    gen_dplus  = max(gen_cdf_diff, [], 2);
    gen_dminus = min(gen_cdf_diff, [], 2);
end

