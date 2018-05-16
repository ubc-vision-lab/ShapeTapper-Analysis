function [obs_out, gen_out] = parse_cdf_struct(cdf_struct)
%PARSECDFSTRUCT Summary of this function goes here
%   Detailed explanation goes here
    obs_out = struct();
    gen_out = struct();
    
    unif_cdf = cdf_struct{1,1}{1,1};

    obs_cdf  = cdf_struct{1,1}{1,2};
    obs_cdf_obs_diff = obs_cdf - unif_cdf;
    [obs_dplus, obs_dplus_r]   = max(obs_cdf_obs_diff);
    [obs_dminus, obs_dminus_r] = min(obs_cdf_obs_diff);
    
    obs_out.dplus   = obs_dplus;
    obs_out.dplus_r = obs_dplus_r;
    obs_out.dminus   = obs_dminus;
    obs_out.dminus_r = obs_dminus_r;
    
    gen_cdf = cdf_struct{1,1}{1,3};
    gen_cdf_diff = gen_cdf - unif_cdf;
    gen_dplus  = max(gen_cdf_diff, [], 2);
    gen_dminus = min(gen_cdf_diff, [], 2);
    
    gen_out.dplus  = gen_dplus;
    gen_out.dplus_avg = mean(gen_dplus);
    gen_out.dplus_std = std(gen_dplus);
    
    gen_out.dminus = gen_dminus;
    gen_out.dminus_avg = mean(gen_dminus);
    gen_out.dminus_std = std(gen_dminus);
end

