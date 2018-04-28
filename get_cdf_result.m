function result = get_cdf_result(dplus, dminus, dplus_r, dminus_r)
%get_cdf_result Returns interpreted point distribution from D+ and D- values
% 'N'=near, 'F'=far, 'Ø'=null(uniform), 'NF'=some near,some far, 'M'=middle
    if dplus > 0.95
            if dminus >= 0.05
                result = {'N'};
            elseif dminus < 0.05 
                if dplus_r <= dminus_r
                    result = {'NF'};
                else 
                    result = {'M'};
                end
            end
    elseif dminus < 0.05
        result = {'F'};
    else 
        result = {'Ø'};
    end
end

