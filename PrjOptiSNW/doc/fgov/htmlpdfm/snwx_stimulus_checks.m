%% Existing Stimulus as a Function of Income and Family Status
% Taking advantage of <https://github.com/FanWangEcon/PrjOptiSNW/blob/master/PrjOptiSNW/fgov/snw_stimulus_checks.m 
% *snw_stimulus_checks*> from the <https://fanwangecon.github.io/PrjOptiSNW/ *PrjOptiSNW 
% Package*,> this function presents stimulus checks at different income levels 
% for households with different children count and martial status. The function 
% considers the first as well as the second stimulus check.
%% Trump Stimulus Checks for Unmarried Households
% Check base amount per adult and per child for the first and second rounds. 

[fl_stimulus_adult_first, fl_stimulus_child_first] = deal(1200, 500);
[fl_stimulus_adult_second, fl_stimulus_child_second] = deal(600, 600);
bl_visualize = true;
%% 
% Visualize stimulus check amounts.

bl_marital = 0;
for it_kids=0:1:4
    snw_stimulus_checks(it_kids, bl_marital, ...
        fl_stimulus_adult_first, fl_stimulus_child_first, ...
        fl_stimulus_adult_second, fl_stimulus_child_second, ...
        bl_visualize);
end
%% Trump Stimulus Checks for arried Households
% Visualize stimulus check amounts.

bl_marital = 1;
for it_kids=0:1:4
    snw_stimulus_checks(it_kids, bl_marital, ...
        fl_stimulus_adult_first, fl_stimulus_child_first, ...
        fl_stimulus_adult_second, fl_stimulus_child_second, ...
        bl_visualize);
end
%% Biden Stimulus Checks for Unmarried and Married Households
% Biden check, what is the maximum phase out given 4 kids and married? 

(1400*6)/(5/100)+150000
%% 
% Check levels (no second round).

[fl_stimulus_adult_first, fl_stimulus_child_first] = deal(1400, 1400);
[fl_stimulus_adult_second, fl_stimulus_child_second] = deal(0, 0);
bl_visualize = true;
%% 
% Visualize stimulus check amounts.

for bl_marital=0:1
    for it_kids=0:1:4
        snw_stimulus_checks(it_kids, bl_marital, ...
            fl_stimulus_adult_first, fl_stimulus_child_first, ...
            fl_stimulus_adult_second, fl_stimulus_child_second, ...
            bl_visualize);
    end
end