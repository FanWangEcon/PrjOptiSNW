%% SNW_STIMULUS_CHECKS_BIDEN Biden checks, not used in Matlab
%    Biden checks, not ujsed in matlab
%
%    See also SNWX_STIMULUS_CHECKS
%

function [varargout] = snw_stimulus_checks_biden(varargin)
% Load inputs or not
if (~isempty(varargin))
    % this relies on externally generated parameters, defaults do not have to be generated
    % if this file has to be invoked many times, then this saves time by avoiding
    % regenerating defaults over and over again
    
    if (length(varargin) ==3)
        [it_kids, bl_marital, bl_visualize] = varargin{:};
    end
    
else
    close all;
    
    bl_visualize = true;
    [it_kids, bl_marital] = deal(0,0);
    
end

% Second round check
fc_stimulus_biden = snwi_stimulus_check(it_kids, bl_marital);

% Visualize
if (bl_visualize)
    figure();
    hold on;
    
    % Set bounds on the domain
    fl_x_lower = 0;
    fl_x_higher = 220000;
    
    % Graph
    pl_check = fplot(fc_stimulus_biden, [fl_x_lower, fl_x_higher],...
        '-', 'Linewidth',2);
    
    % Legend
    legend([pl_check], ...
        {'Biden checks'}, ...
        'Location','best',...
        'NumColumns',1,'FontSize',12,'TextColor','black');
    
    % Add x-axis and y-axis
    pl_xline = xline(0);
    pl_yline = yline(0);
    pl_xline.HandleVisibility = 'off';
    pl_yline.HandleVisibility = 'off';
    
    % Title and y and y-able
    title({'Biden stimulus check amounts as a function of household income (under consideration)',...
        ['Kids count = ' num2str(it_kids) '; marital status = ' num2str(bl_marital)]}, ...
        'Interpreter', "none");
    ylabel('Stimulus check amount');
    xlabel('income levels');
    
    % Add grids
    grid on;
    grid minor;
end

%% Return 
varargout = cell(nargout,0);
for it_k = 1:nargout
    if (it_k==1)
        ob_out_cur = fc_stimulus_biden;
    end
    varargout{it_k} = ob_out_cur;
end

end

%% Biden Check Actual
function fc_stimulus_check = snwi_stimulus_check(it_kids, bl_marital)

fl_stimulus_adult = 1400;
fl_stimulus_child = 1400;
fl_slope_m0k0 = 1400/(80000-75000);
fl_slope_m0k1 = 2800/(120000-112500);
fl_slope_m0k2 = 4200/(120000-112500);
fl_slope_m0k3 = 5600/(120000-112500);
fl_slope_m0k4 = 7000/(120000-112500);
% disp(['Check Reduction Per $100 fl_slope_m0k0=' num2str(round(fl_slope_m0k0*100, 3))]);
% disp(['Check Reduction Per $100 fl_slope_m0k1=' num2str(round(fl_slope_m0k1*100, 3))]);
% disp(['Check Reduction Per $100 fl_slope_m0k2=' num2str(round(fl_slope_m0k2*100, 3))]);
% disp(['Check Reduction Per $100 fl_slope_m0k3=' num2str(round(fl_slope_m0k3*100, 3))]);
% disp(['Check Reduction Per $100 fl_slope_m0k4=' num2str(round(fl_slope_m0k4*100, 3))]);

fl_slope_m1k0 = 2800/(160000-150000);
fl_slope_m1k1 = 4200/(160000-150000);
fl_slope_m1k2 = 5600/(160000-150000);
fl_slope_m1k3 = 7000/(160000-150000);
fl_slope_m1k4 = 8400/(160000-150000);
% disp(['Check Reduction Per $100 fl_slope_m1k0=' num2str(round(fl_slope_m1k0*100, 3))]);
% disp(['Check Reduction Per $100 fl_slope_m1k1=' num2str(round(fl_slope_m1k1*100, 3))]);
% disp(['Check Reduction Per $100 fl_slope_m1k2=' num2str(round(fl_slope_m1k2*100, 3))]);
% disp(['Check Reduction Per $100 fl_slope_m1k3=' num2str(round(fl_slope_m1k3*100, 3))]);
% disp(['Check Reduction Per $100 fl_slope_m1k4=' num2str(round(fl_slope_m1k4*100, 3))]);

% Check Reduction Per $100 fl_slope_m0k0=28
% Check Reduction Per $100 fl_slope_m0k1=37.333
% Check Reduction Per $100 fl_slope_m0k2=56
% Check Reduction Per $100 fl_slope_m0k3=74.667
% Check Reduction Per $100 fl_slope_m0k4=93.333
% Check Reduction Per $100 fl_slope_m1k0=28
% Check Reduction Per $100 fl_slope_m1k1=42
% Check Reduction Per $100 fl_slope_m1k2=56
% Check Reduction Per $100 fl_slope_m1k3=70
% Check Reduction Per $100 fl_slope_m1k4=84

% start point household head
fl_check_dollar = fl_stimulus_adult;
% married household gets more, marital = 0, 1
fl_check_dollar = fl_check_dollar + bl_marital*fl_stimulus_adult;
% Households with kids: 0, 1,2,3,4
fl_check_dollar = fl_check_dollar + it_kids*fl_stimulus_child;

if (~bl_marital && it_kids == 0)
    % The benefit would start decreasing at a rate of $5 for every additional $100 in income
    fc_check_reduce = @(x) ((max(x - 75000,0))*fl_slope_m0k0);
end

% phaseout starts $112,500 for heads of household
if (~bl_marital && it_kids == 1)    
    fc_check_reduce = @(x) ((max(x - 112500,0))*fl_slope_m0k1);
elseif (~bl_marital && it_kids == 2)
    fc_check_reduce = @(x) ((max(x - 112500,0))*fl_slope_m0k2);
elseif (~bl_marital && it_kids == 3)
    fc_check_reduce = @(x) ((max(x - 112500,0))*fl_slope_m0k3);
elseif (~bl_marital && it_kids == 4)
    fc_check_reduce = @(x) ((max(x - 112500,0))*fl_slope_m0k4);
end

% phaseout starts $150,000 for heads of household
if (bl_marital && it_kids == 0)
    fc_check_reduce = @(x) ((max(x - 150000,0))*fl_slope_m1k0);
elseif (bl_marital && it_kids == 1)
    fc_check_reduce = @(x) ((max(x - 150000,0))*fl_slope_m1k1);
elseif (bl_marital && it_kids == 2)
    fc_check_reduce = @(x) ((max(x - 150000,0))*fl_slope_m1k2);
elseif (bl_marital && it_kids == 3)
    fc_check_reduce = @(x) ((max(x - 150000,0))*fl_slope_m1k3);
elseif (bl_marital && it_kids == 4)
    fc_check_reduce = @(x) ((max(x - 150000,0))*fl_slope_m1k4);
end

% stimulus check function first check
fc_stimulus_check = @(x) max(0, fl_check_dollar - fc_check_reduce(x));

end

%% Biden Check Proposed
function fc_stimulus_check = snwi_stimulus_check_proposed(it_kids, bl_marital)
% x is household income
% ,Single with 0 children,Stimulus check,,,Single with 1 child,Stimulus check,,,Single with 2 children,Stimulus check,,,Single with 3 children,Stimulus check,,,Single with 4 children,Stimulus check,,,,,Married,Child count
% Income,75000,1400,,Income,112500,2800,,Income,112500,4200,,Income,112500,5600,,Income,112500,7000,,,,,,
% Income,100000,0,,Income,150000,0,,Income,150000,0,,Income,150000,0,,Income,150000,0,,,,,,
% ,Married with 0 children,Stimulus check,,,Married with 1 child,Stimulus check,,,Married with 2 children,Stimulus check,,,Married with 3 children,Stimulus check,,,Married with 4 children,Stimulus check,,,,,,
% Income,150000,2800,,Income,150000,4200,,Income,150000,5600,,Income,150000,7000,,Income,150000,8400,,,,,,
% Income,200000,0,,Income,200000,0,,Income,200000,0,,Income,200000,0,,Income,200000,0,,,,,,

fl_stimulus_adult = 1400;
fl_stimulus_child = 1400;
fl_slope_m0k0 = 1400/(100000-75000);
fl_slope_m0k1 = 2800/(150000-112500);
fl_slope_m0k2 = 4200/(150000-112500);
fl_slope_m0k3 = 5600/(150000-112500);
fl_slope_m0k4 = 7000/(150000-112500);
% disp(['Check Reduction Per $100 fl_slope_m0k0=' num2str(round(fl_slope_m0k0*100, 3))]);
% disp(['Check Reduction Per $100 fl_slope_m0k1=' num2str(round(fl_slope_m0k1*100, 3))]);
% disp(['Check Reduction Per $100 fl_slope_m0k2=' num2str(round(fl_slope_m0k2*100, 3))]);
% disp(['Check Reduction Per $100 fl_slope_m0k3=' num2str(round(fl_slope_m0k3*100, 3))]);
% disp(['Check Reduction Per $100 fl_slope_m0k4=' num2str(round(fl_slope_m0k4*100, 3))]);

fl_slope_m1k0 = 2800/(200000-150000);
fl_slope_m1k1 = 4200/(200000-150000);
fl_slope_m1k2 = 5600/(200000-150000);
fl_slope_m1k3 = 7000/(200000-150000);
fl_slope_m1k4 = 8400/(200000-150000);
% disp(['Check Reduction Per $100 fl_slope_m1k0=' num2str(round(fl_slope_m1k0*100, 3))]);
% disp(['Check Reduction Per $100 fl_slope_m1k1=' num2str(round(fl_slope_m1k1*100, 3))]);
% disp(['Check Reduction Per $100 fl_slope_m1k2=' num2str(round(fl_slope_m1k2*100, 3))]);
% disp(['Check Reduction Per $100 fl_slope_m1k3=' num2str(round(fl_slope_m1k3*100, 3))]);
% disp(['Check Reduction Per $100 fl_slope_m1k4=' num2str(round(fl_slope_m1k4*100, 3))]);

% Check Reduction Per $100 fl_slope_m0k0=5.6
% Check Reduction Per $100 fl_slope_m0k1=7.467
% Check Reduction Per $100 fl_slope_m0k2=11.2
% Check Reduction Per $100 fl_slope_m0k3=14.933
% Check Reduction Per $100 fl_slope_m0k4=18.667
% Check Reduction Per $100 fl_slope_m1k0=5.6
% Check Reduction Per $100 fl_slope_m1k1=8.4
% Check Reduction Per $100 fl_slope_m1k2=11.2
% Check Reduction Per $100 fl_slope_m1k3=14
% Check Reduction Per $100 fl_slope_m1k4=16.8

% start point household head
fl_check_dollar = fl_stimulus_adult;
% married household gets more, marital = 0, 1
fl_check_dollar = fl_check_dollar + bl_marital*fl_stimulus_adult;
% Households with kids: 0, 1,2,3,4
fl_check_dollar = fl_check_dollar + it_kids*fl_stimulus_child;

if (~bl_marital && it_kids == 0)
    % The benefit would start decreasing at a rate of $5 for every additional $100 in income
    fc_check_reduce = @(x) ((max(x - 75000,0))*fl_slope_m0k0);
end

% phaseout starts $112,500 for heads of household
if (~bl_marital && it_kids == 1)    
    fc_check_reduce = @(x) ((max(x - 112500,0))*fl_slope_m0k1);
elseif (~bl_marital && it_kids == 2)
    fc_check_reduce = @(x) ((max(x - 112500,0))*fl_slope_m0k2);
elseif (~bl_marital && it_kids == 3)
    fc_check_reduce = @(x) ((max(x - 112500,0))*fl_slope_m0k3);
elseif (~bl_marital && it_kids == 4)
    fc_check_reduce = @(x) ((max(x - 112500,0))*fl_slope_m0k4);
end

% phaseout starts $150,000 for heads of household
if (bl_marital && it_kids == 0)
    fc_check_reduce = @(x) ((max(x - 150000,0))*fl_slope_m1k0);
elseif (bl_marital && it_kids == 1)
    fc_check_reduce = @(x) ((max(x - 150000,0))*fl_slope_m1k1);
elseif (bl_marital && it_kids == 2)
    fc_check_reduce = @(x) ((max(x - 150000,0))*fl_slope_m1k2);
elseif (bl_marital && it_kids == 3)
    fc_check_reduce = @(x) ((max(x - 150000,0))*fl_slope_m1k3);
elseif (bl_marital && it_kids == 4)
    fc_check_reduce = @(x) ((max(x - 150000,0))*fl_slope_m1k4);
end

% stimulus check function first check
fc_stimulus_check = @(x) max(0, fl_check_dollar - fc_check_reduce(x));

end
