%% Existing Stimulus as a Function of Income and Family Status
% Taking advantage of <https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/fgov/snw_stimulus_checks_biden.m 
% *snw_stimulus_checks_biden*> from the <https://fanwangecon.github.io/PrjOptiSNW/ 
% *PrjOptiSNW Package*,> this function presents stimulus checks at different income 
% levels for households with different children count and martial status. This 
% is for the Biden Stimulus under the <https://en.wikipedia.org/wiki/American_Rescue_Plan_Act_of_2021#:~:text=No.,pandemic%20and%20the%20ongoing%20recession. 
% American Rescue Plan Act>.
%% Biden Stimulus Checks for Unmarried Households
% Check base amount per adult and per child for the first and second rounds. 
% 
% Visualize stimulus check amounts.

bl_visualize = true;
bl_marital = 0;
for it_kids=0:1:4
    snw_stimulus_checks_biden(it_kids, bl_marital, bl_visualize);
end
%% Biden Stimulus Checks for Married Households
% Visualize stimulus check amounts.

bl_marital = 1;
for it_kids=0:1:4
    snw_stimulus_checks_biden(it_kids, bl_marital, bl_visualize);
end
%%