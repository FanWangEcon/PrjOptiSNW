function welfare = snw_hh_welfare(V, gamma)
% Convert from lifetime V to Welfare in consumption units
% V convertion is based on finding constant life-time consumption flow that
% will equal to V, then doing life-discounting based adjustments.

welfare = (V.*(1-gamma)).^(1/(1-gamma));
welfare(V==0) = 0;

end
