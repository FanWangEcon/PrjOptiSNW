function F=spousal_income(j,educ,kids,earn,SS_inc)

global jret

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
        F=coeff_cons+(coeff_ln_Rincome*log(earn+SS_inc))+(coeff_Rage*(j+17))+(coeff_Rage_sq*((j+17)^2))+(coeff_Rage_cube*((j+17)^3))+(coeff_Rcollege*(educ-1));
    elseif kids==2
        F=coeff_cons+(coeff_ln_Rincome*log(earn+SS_inc))+(coeff_Rage*(j+17))+(coeff_Rage_sq*((j+17)^2))+(coeff_Rage_cube*((j+17)^3))+(coeff_Rcollege*(educ-1))+coeff_Rkids1;
    elseif kids==3
        F=coeff_cons+(coeff_ln_Rincome*log(earn+SS_inc))+(coeff_Rage*(j+17))+(coeff_Rage_sq*((j+17)^2))+(coeff_Rage_cube*((j+17)^3))+(coeff_Rcollege*(educ-1))+coeff_Rkids2;
    elseif kids==4
        F=coeff_cons+(coeff_ln_Rincome*log(earn+SS_inc))+(coeff_Rage*(j+17))+(coeff_Rage_sq*((j+17)^2))+(coeff_Rage_cube*((j+17)^3))+(coeff_Rcollege*(educ-1))+coeff_Rkids3;
    elseif kids==5
        F=coeff_cons+(coeff_ln_Rincome*log(earn+SS_inc))+(coeff_Rage*(j+17))+(coeff_Rage_sq*((j+17)^2))+(coeff_Rage_cube*((j+17)^3))+(coeff_Rcollege*(educ-1))+coeff_Rkids4;
    else
        error('Update spousal_income.m to allow for more than 4 kids')
    end
else
    if kids==1
        F=coeff_cons+(coeff_ln_Rincome*log(earn+SS_inc))+(coeff_Rage*(j+17))+(coeff_Rage_sq*((j+17)^2))+(coeff_Rage_cube*((j+17)^3))+(coeff_Rcollege*(educ-1))+coeff_dummy65+(coeff_dummy65kids*(kids-1));
    elseif kids==2
        F=coeff_cons+(coeff_ln_Rincome*log(earn+SS_inc))+(coeff_Rage*(j+17))+(coeff_Rage_sq*((j+17)^2))+(coeff_Rage_cube*((j+17)^3))+(coeff_Rcollege*(educ-1))+coeff_dummy65+(coeff_dummy65kids*(kids-1))+coeff_Rkids1;
    elseif kids==3
        F=coeff_cons+(coeff_ln_Rincome*log(earn+SS_inc))+(coeff_Rage*(j+17))+(coeff_Rage_sq*((j+17)^2))+(coeff_Rage_cube*((j+17)^3))+(coeff_Rcollege*(educ-1))+coeff_dummy65+(coeff_dummy65kids*(kids-1))+coeff_Rkids2;
    elseif kids==4
        F=coeff_cons+(coeff_ln_Rincome*log(earn+SS_inc))+(coeff_Rage*(j+17))+(coeff_Rage_sq*((j+17)^2))+(coeff_Rage_cube*((j+17)^3))+(coeff_Rcollege*(educ-1))+coeff_dummy65+(coeff_dummy65kids*(kids-1))+coeff_Rkids3;
    elseif kids==5
        F=coeff_cons+(coeff_ln_Rincome*log(earn+SS_inc))+(coeff_Rage*(j+17))+(coeff_Rage_sq*((j+17)^2))+(coeff_Rage_cube*((j+17)^3))+(coeff_Rcollege*(educ-1))+coeff_dummy65+(coeff_dummy65kids*(kids-1))+coeff_Rkids4;
    else
        error('Update spousal_income.m to allow for more than 4 kids')
    end
end

F=exp(F);

% if j<jret
%
%     % reg ln_Sincome ln_Rincome c.RAGE##c.RAGE##c.RAGE Rcollege i.kids dummy_YEAR* if inrange(RAGE,18,64) & Rmarried==1 [pweight=WEIGHT]
%     % Income refers to labor income+Social Security+SSI+UI benefits+other transfers
%     coeff_ln_Rincome=0.1110694;
%     coeff_Rage=0.251775;
%     coeff_Rage_sq=-0.004877;
%     coeff_Rage_cube=0.0000291;
%     coeff_Rcollege=0.2245369;
%     coeff_Rkids1=-0.153976;
%     coeff_Rkids2=-0.3002043;
%     coeff_Rkids3=-0.4702531;
%     coeff_Rkids4=-0.7895856;
%     coeff_cons=-4.58432;
%
%     if kids==1
%         F=coeff_cons+(coeff_ln_Rincome*log(earn+SS_inc))+(coeff_Rage*(j+17))+(coeff_Rage_sq*((j+17)^2))+(coeff_Rage_cube*((j+17)^3))+(coeff_Rcollege*(educ-1));
%     elseif kids==2
%         F=coeff_cons+(coeff_ln_Rincome*log(earn+SS_inc))+(coeff_Rage*(j+17))+(coeff_Rage_sq*((j+17)^2))+(coeff_Rage_cube*((j+17)^3))+(coeff_Rcollege*(educ-1))+coeff_Rkids1;
%     elseif kids==3
%         F=coeff_cons+(coeff_ln_Rincome*log(earn+SS_inc))+(coeff_Rage*(j+17))+(coeff_Rage_sq*((j+17)^2))+(coeff_Rage_cube*((j+17)^3))+(coeff_Rcollege*(educ-1))+coeff_Rkids2;
%     elseif kids==4
%         F=coeff_cons+(coeff_ln_Rincome*log(earn+SS_inc))+(coeff_Rage*(j+17))+(coeff_Rage_sq*((j+17)^2))+(coeff_Rage_cube*((j+17)^3))+(coeff_Rcollege*(educ-1))+coeff_Rkids3;
%     elseif kids==5
%         F=coeff_cons+(coeff_ln_Rincome*log(earn+SS_inc))+(coeff_Rage*(j+17))+(coeff_Rage_sq*((j+17)^2))+(coeff_Rage_cube*((j+17)^3))+(coeff_Rcollege*(educ-1))+coeff_Rkids4;
%     else
%         error('Update spousal_income.m to allow for more than 4 kids')
%     end
% else
%
%     % reg ln_Sincome ln_Rincome c.RAGE Rcollege dummy_YEAR* if RAGE>=65 & RAGE!=. & Rmarried==1 [pweight=WEIGHT]
%     % Income refers to labor income+Social Security+SSI+UI benefits+other transfers
%     coeff_ln_Rincome=0.2213816;
%     coeff_Rage=-0.0305283;
%     coeff_Rcollege=0.1391068;
%     coeff_cons=0.7526077;
%
%     F=coeff_cons+(coeff_ln_Rincome*log(earn+SS_inc))+(coeff_Rage*(j+17))+(coeff_Rcollege*(educ-1));
% end
%
% F=exp(F);

end
