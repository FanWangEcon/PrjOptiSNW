%% SNW_STIMULUS_CHECKS provides the amount of stimulus checks per household type
%    Considers household child count, martial status. Function provides
%    stimulus check amount in dollars given household income level. Using
%    US dollar as unit of accounting. Provides the first check under Trump,
%    the second check under Trump, and both checks under Trump jointly. 
%
%    * see: <https://www.irs.gov/pub/irs-utl/e-poster_payments.pdf first
%    check>
%    * see:
%    <https://www.irs.gov/newsroom/treasury-and-irs-begin-delivering-second-round-of-economic-impact-payments-to-millions-of-americans,
%    second check>
%    * see:
%    <https://github.com/FanWangEcon/PrjOptiAlloc/blob/master/R/ffp_snw_welfarechecks_input.R
%    ffp_snw_welfarechecks_input>
%
%    See also SNWX_STIMULUS_CHECKS
%

function [varargout] = snw_stimulus_checks(varargin)

% Defaults
[fl_stimulus_adult_first, fl_stimulus_child_first] = deal(1200, 500);
[fl_stimulus_adult_second, fl_stimulus_child_second] = deal(600, 600);

% Load inputs or not
if (~isempty(varargin))
    % this relies on externally generated parameters, defaults do not have to be generated
    % if this file has to be invoked many times, then this saves time by avoiding
    % regenerating defaults over and over again
    
    bl_visualize = false;
    if (length(varargin) ==2)
        [it_kids, bl_marital] = varargin{:};
    elseif (length(varargin) == 6)
        [it_kids, bl_marital, ...
            fl_stimulus_adult_first, fl_stimulus_child_first, ...
            fl_stimulus_adult_second, fl_stimulus_child_second] = varargin{:};
    elseif (length(varargin) == 7)
        [it_kids, bl_marital, ...
            fl_stimulus_adult_first, fl_stimulus_child_first, ...
            fl_stimulus_adult_second, fl_stimulus_child_second, ...
            bl_visualize] = varargin{:};
    end
    
else
    close all;
    
    bl_visualize = true;
    [it_kids, bl_marital] = deal(0,0);
    
end

% First round check
fc_stimulus_check_first = snwi_stimulus_check_trump(it_kids, bl_marital, ...
    fl_stimulus_adult_first, fl_stimulus_child_first);

% Second round check
fc_stimulus_check_second = snwi_stimulus_check_trump(it_kids, bl_marital, ...
    fl_stimulus_adult_second, fl_stimulus_child_second);

% Both checks together
fc_stimulus_check_firstsecond = @(x) fc_stimulus_check_first(x) + fc_stimulus_check_second(x);

% Visualize
if (bl_visualize)
    figure();
    hold on;
    
    % Set bounds on the domain
    fl_x_lower = 0;
    fl_x_higher = 350000;
    
    % Graph
    pl_check1 = fplot(fc_stimulus_check_first, [fl_x_lower, fl_x_higher], ...
        '-', 'Linewidth',10);
    pl_check2 = fplot(fc_stimulus_check_second, [fl_x_lower, fl_x_higher], ...
        '-', 'Linewidth',5);
    pl_bothchecks = fplot(fc_stimulus_check_firstsecond, [fl_x_lower, fl_x_higher],...
        '-', 'Linewidth',2);
    
    % Legend
    legend([pl_check1, pl_check2, pl_bothchecks], ...
        {'First check','Second check','Both checks'}, ...
        'Location','best',...
        'NumColumns',1,'FontSize',12,'TextColor','black');
    
    % Add x-axis and y-axis
    pl_xline = xline(0);
    pl_yline = yline(0);
    pl_xline.HandleVisibility = 'off';
    pl_yline.HandleVisibility = 'off';
    
    % Title and y and y-able
    title({'Stimulus Check amounts as a function of household income',...
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
        ob_out_cur = fc_stimulus_check_firstsecond;
    elseif (it_k==2)
        ob_out_cur = fc_stimulus_check_first;
    elseif (it_k==3)
        ob_out_cur = fc_stimulus_check_second;
    end
    varargout{it_k} = ob_out_cur;
end

end

%% Trump Check
function fc_stimulus_check = snwi_stimulus_check_trump(...
    it_kids, bl_marital, fl_stimulus_adult, fl_stimulus_child)
% x is household income
% fl_stimulus_adult = 1200, fl_stimulus_adult = 500 for first round
% fl_stimulus_adult = 600, fl_stimulus_adult = 600 for second round

% start point household head
fl_check_dollar = fl_stimulus_adult;
% married household gets more, marital = 0, 1
fl_check_dollar = fl_check_dollar + bl_marital*fl_stimulus_adult;
% Households with kids: 0, 1,2,3,4
fl_check_dollar = fl_check_dollar + it_kids*fl_stimulus_child;

if (~bl_marital && it_kids == 0)
    % The benefit would start decreasing at a rate of $5 for every additional $100 in income
    fc_check_reduce = @(x) ((max(x - 75000,0))/100)*5;
end

% phaseout starts $112,500 for heads of household
if (~bl_marital && it_kids ~= 0)
    % The benefit would start decreasing at a rate of $5 for every additional $100 in income
    fc_check_reduce = @(x) ((max(x - 112500,0))/100)*5;
end

% phaseout starts $150,000 for heads of household
if (bl_marital)
    % The benefit would start decreasing at a rate of $5 for every additional $100 in income
    fc_check_reduce = @(x) ((max(x - 150000,0))/100)*5;
end

% stimulus check function first check
fc_stimulus_check = @(x) max(0, fl_check_dollar - fc_check_reduce(x));

end