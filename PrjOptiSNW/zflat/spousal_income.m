function F=spousal_income(j,educ,kids,earn,SS_inc)

global jret
 
% if j<jret
%     
%     % reg ln_SlabSS ln_RLABINC RAGE RAGE_sq RAGE_cube Rcollege kids dummy_YEAR1 dummy_YEAR2 dummy_YEAR3 dummy_YEAR4 dummy_YEAR5 dummy_YEAR6 dummy_YEAR7 dummy_YEAR8 dummy_YEAR9 dummy_YEAR10 dummy_YEAR11 if inrange(RAGE,18,65) & Rmarried==1 [pweight=WEIGHT]
%     coeff_Rearn=0.1141646; % 0.1021238; %0.1000593;
%     coeff_Rage=0.2366304; % 0.2577159; % 0.2769428;
%     coeff_Rage_sq=-0.0044643; % -0.0049819; % -0.0054628;
%     coeff_Rage_cube=0.0000255; % 0.0000296; % 0.0000337;
%     coeff_Reduc=0.2241748; % 0.2231719; % 0.2140897;
%     coeff_Rkids=-0.1670697; % -0.1674179; % -0.1685701; % set equal to 0 to take out kid effect
%     coeff_cons=-4.606553; % -5.039279; % -5.035496;
% 
%     F=coeff_cons+(coeff_Rearn*log(earn))+(coeff_Rage*(j+17))+(coeff_Rage_sq*((j+17)^2))+(coeff_Rage_cube*((j+17)^3))+(coeff_Reduc*(educ-1))+(coeff_Rkids*(kids-1));
% 
%     F=exp(F);
%     
% else
%     
%     %  reg ln_SlabSS ln_RlabSS RAGE Rcollege dummy_YEAR1 dummy_YEAR2 dummy_YEAR3 dummy_YEAR4 dummy_YEAR5 dummy_YEAR6 dummy_YEAR7 dummy_YEAR8 dummy_YEAR9 dummy_YEAR10 dummy_YEAR11 if RAGE>=66 & RAGE!=. & Rmarried==1 [pweight=WEIGHT]
%     coeff_RSSinc=0.2275314;  % 0.2929644;
%     coeff_Rage=-0.0291802;  % -0.0345724;
%     coeff_Reduc=0.1444211;  % 0.1771493;
%     coeff_cons=0.537291;  % -2.012454;
% 
%     F=coeff_cons+(coeff_RSSinc*log(SS_inc))+(coeff_Rage*(j+17))+(coeff_Reduc*(educ-1));
% 
%     F=exp(F);
%     
% end

% reg ln_Sincome ln_Rincome c.RAGE##c.RAGE##c.RAGE Rcollege i.kids dummy65 dummy65kids dummy_YEAR* if RAGE>=18 & RAGE!=. & Rmarried==1 [pweight=WEIGHT]
% Income refers to labor income+Social Security+SSI+UI benefits+other transfers
coeff_ln_Rincome=0.1344672;
coeff_Rage=0.1639199;
coeff_Rage_sq=-0.0027544;
coeff_Rage_cube=0.0000127;
coeff_Rcollege=0.2061415;
coeff_Rkids1=-0.1675939;
coeff_Rkids2=-0.3018777;
coeff_Rkids3=-0.4624552;
coeff_Rkids4=-0.7861216;
coeff_dummy65=-0.222999;
coeff_dummy65kids=0.1727698;
coeff_cons=-3.407486;

if j<jret
    if kids==1
        F=coeff_cons+(coeff_ln_Rincome*log(earn))+(coeff_Rage*(j+17))+(coeff_Rage_sq*((j+17)^2))+(coeff_Rage_cube*((j+17)^3))+(coeff_Rcollege*(educ-1));
    elseif kids==2
        F=coeff_cons+(coeff_ln_Rincome*log(earn))+(coeff_Rage*(j+17))+(coeff_Rage_sq*((j+17)^2))+(coeff_Rage_cube*((j+17)^3))+(coeff_Rcollege*(educ-1))+coeff_Rkids1;
    elseif kids==3
        F=coeff_cons+(coeff_ln_Rincome*log(earn))+(coeff_Rage*(j+17))+(coeff_Rage_sq*((j+17)^2))+(coeff_Rage_cube*((j+17)^3))+(coeff_Rcollege*(educ-1))+coeff_Rkids2;
    elseif kids==4
        F=coeff_cons+(coeff_ln_Rincome*log(earn))+(coeff_Rage*(j+17))+(coeff_Rage_sq*((j+17)^2))+(coeff_Rage_cube*((j+17)^3))+(coeff_Rcollege*(educ-1))+coeff_Rkids3;
    elseif kids==5
        F=coeff_cons+(coeff_ln_Rincome*log(earn))+(coeff_Rage*(j+17))+(coeff_Rage_sq*((j+17)^2))+(coeff_Rage_cube*((j+17)^3))+(coeff_Rcollege*(educ-1))+coeff_Rkids4;
    else
        error('Update spousal_income.m to allow for more than 4 kids')
    end
else
    if kids==1
        F=coeff_cons+(coeff_ln_Rincome*log(earn))+(coeff_Rage*(j+17))+(coeff_Rage_sq*((j+17)^2))+(coeff_Rage_cube*((j+17)^3))+(coeff_Rcollege*(educ-1))+coeff_dummy65+(coeff_dummy65kids*(kids-1));
    elseif kids==2
        F=coeff_cons+(coeff_ln_Rincome*log(earn))+(coeff_Rage*(j+17))+(coeff_Rage_sq*((j+17)^2))+(coeff_Rage_cube*((j+17)^3))+(coeff_Rcollege*(educ-1))+coeff_dummy65+(coeff_dummy65kids*(kids-1))+coeff_Rkids1;
    elseif kids==3
        F=coeff_cons+(coeff_ln_Rincome*log(earn))+(coeff_Rage*(j+17))+(coeff_Rage_sq*((j+17)^2))+(coeff_Rage_cube*((j+17)^3))+(coeff_Rcollege*(educ-1))+coeff_dummy65+(coeff_dummy65kids*(kids-1))+coeff_Rkids2;
    elseif kids==4
        F=coeff_cons+(coeff_ln_Rincome*log(earn))+(coeff_Rage*(j+17))+(coeff_Rage_sq*((j+17)^2))+(coeff_Rage_cube*((j+17)^3))+(coeff_Rcollege*(educ-1))+coeff_dummy65+(coeff_dummy65kids*(kids-1))+coeff_Rkids3;
    elseif kids==5
        F=coeff_cons+(coeff_ln_Rincome*log(earn))+(coeff_Rage*(j+17))+(coeff_Rage_sq*((j+17)^2))+(coeff_Rage_cube*((j+17)^3))+(coeff_Rcollege*(educ-1))+coeff_dummy65+(coeff_dummy65kids*(kids-1))+coeff_Rkids4;
    else
        error('Update spousal_income.m to allow for more than 4 kids')
    end
end

F=exp(F);

end