function F=spousal_income(j,educ,kids,earn,SS_inc)

global jret

if j<jret
    
    % reg ln_SlabSS ln_RLABINC RAGE RAGE_sq RAGE_cube Rcollege kids dummy_YEAR1 dummy_YEAR2 dummy_YEAR3 dummy_YEAR4 dummy_YEAR5 dummy_YEAR6 dummy_YEAR7 dummy_YEAR8 dummy_YEAR9 dummy_YEAR10 dummy_YEAR11 if inrange(RAGE,18,65) & Rmarried==1 [pweight=WEIGHT]
    coeff_Rearn=0.1021238; %0.1000593;
    coeff_Rage=0.2577159; % 0.2769428;
    coeff_Rage_sq=-0.0049819; % -0.0054628;
    coeff_Rage_cube=0.0000296; % 0.0000337;
    coeff_Reduc=0.2231719; % 0.2140897;
    coeff_Rkids=-0.1674179; % -0.1685701; % set equal to 0 to take out kid effect
    coeff_cons=-5.039279; % -5.035496;

    F=coeff_cons+(coeff_Rearn*log(earn))+(coeff_Rage*(j+17))+(coeff_Rage_sq*((j+17)^2))+(coeff_Rage_cube*((j+17)^3))+(coeff_Reduc*(educ-1))+(coeff_Rkids*(kids-1));

    F=exp(F);
    
else
    
    %  reg ln_SlabSS ln_RSSINC RAGE Rcollege dummy_YEAR1 dummy_YEAR2 dummy_YEAR3 dummy_YEAR4 dummy_YEAR5 dummy_YEAR6 dummy_YEAR7 dummy_YEAR8 dummy_YEAR9 dummy_YEAR10 dummy_YEAR11 if RAGE>=66 & RAGE!=. & Rmarried==1 [pweight=WEIGHT]
    coeff_RSSinc=0.2929644;
    coeff_Rage=-0.0345724;
    coeff_Reduc=0.1771493;
    coeff_cons=-2.012454;

    F=coeff_cons+(coeff_RSSinc*log(SS_inc))+(coeff_Rage*(j+17))+(coeff_Reduc*(educ-1));

    F=exp(F);
    
end

end