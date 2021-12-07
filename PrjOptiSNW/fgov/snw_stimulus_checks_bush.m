%% SNW_STIMULUS_CHECKS_BUSH Bush checks, demonstrative, not used
%    Bush checks, not used in matlab, meaning matlab solution functions do
%    not call this. There are separate R functions that implement the
%    checks that are unrelated to this, this is just to graph out the
%    stimulus schedule for visualization
%
%    Economic Stimulus Act of 2008: 600 ($1,200) for singles (couples)
%    making less than $75,000 ($150,000). $300 per child. Amount phases out
%    at a rate of 5 percent (drops by $50 per $1,000 in income exceeding
%    $75,000 ($150,000)).
%
%    IRS rules:
%    <https://www.irs.gov/newsroom/single-head-of-household-with-children
%    single with children>,
%    <https://www.irs.gov/newsroom/married-with-children married with
%    children>,
%    <https://www.irs.gov/newsroom/married-without-qualifying-children
%    married no children>, and
%    <https://www.irs.gov/newsroom/single-without-qualifying-children
%    single no children>.
%
%    See also SNWX_STIMULUS_CHECKS_BUSH, SNWX_STIMULUS_CHECKS_BIDEN
%

function [varargout] = snw_stimulus_checks_bush(varargin)

% Defaults
[fl_stimulus_adult, fl_stimulus_child] = deal(600, 300);

if (~isempty(varargin))
    % this relies on externally generated parameters, defaults do not have to be generated
    % if this file has to be invoked many times, then this saves time by avoiding
    % regenerating defaults over and over again

    if (length(varargin) ==3)
        [it_kids, bl_marital, bl_visualize] = varargin{:};
    elseif (length(varargin) ==5)
        [it_kids, bl_marital, ...
            fl_stimulus_adult, fl_stimulus_child, ...
            bl_visualize] = varargin{:};
    end

else
    close all;

    bl_visualize = true;
    [it_kids, bl_marital] = deal(0,0);

end

% Bush 2008 stimulus
% Set bounds on the domain
fl_x_lower = 0;
fl_x_higher = 250000;

ar_x = linspace(fl_x_lower, fl_x_higher, 1000);
ar_stimulus_bush = snwi_stimulus_check_bush(it_kids, bl_marital, fl_stimulus_adult, fl_stimulus_child, ar_x);

% Visualize
if (bl_visualize)
    figure();
    hold on;


    % Graph
    pl_check = plot(ar_x, ar_stimulus_bush,...
        '-', 'Linewidth',2);

    % Legend
    legend([pl_check], ...
        {'Bush 2008 Stimulus Checks'}, ...
        'Location','best',...
        'NumColumns',1,'FontSize',12,'TextColor','black');

    % Add x-axis and y-axis
    pl_xline = xline(0);
    pl_yline = yline(0);
    pl_xline.HandleVisibility = 'off';
    pl_yline.HandleVisibility = 'off';

    % Title and y and y-able
    title({'Bush stimulus check amounts as a function of household income (under consideration)',...
        ['Kids count = ' num2str(it_kids) '; marital status = ' num2str(bl_marital)]}, ...
        'Interpreter', "none");
    ylabel('Stimulus check amount (Economic Stimulus Act of 2008)');
    xlabel('Income levels');

    % Add grids
    grid on;
    grid minor;
end

%% Return
varargout = cell(nargout,0);
for it_k = 1:nargout
    if (it_k==1)
        ob_out_cur = fc_stimulus_bush;
    end
    varargout{it_k} = ob_out_cur;
end

end

%% Bush Check Approximate
function fc_stimulus_check = snwi_stimulus_check_bush_approximate(...
    it_kids, bl_marital, fl_stimulus_adult, fl_stimulus_child)
% x is household income
% fl_stimulus_adult = 600, fl_stimulus_adult = 600 for first round

% start point household head
fl_check_dollar = fl_stimulus_adult;
% married household gets more, marital = 0, 1
fl_check_dollar = fl_check_dollar + bl_marital*fl_stimulus_adult;
% Households with kids: 0, 1,2,3,4
fl_check_dollar = fl_check_dollar + it_kids*fl_stimulus_child;

% https://www.irs.gov/newsroom/single-without-qualifying-children
if (~bl_marital && it_kids == 0)
    % The benefit would start decreasing at a rate of $5 for every additional $100 in income
    fc_check_reduce = @(x) ((max(x - 75000,0))/100)*5;
end

% phaseout starts $75000 for heads of household
% https://www.irs.gov/newsroom/single-head-of-household-with-children
if (~bl_marital && it_kids ~= 0)
    % The benefit would start decreasing at a rate of $5 for every additional $100 in income
    fc_check_reduce = @(x) ((max(x - 75000,0))/100)*5;
end

% phaseout starts $150,000 for heads of household
% https://www.irs.gov/newsroom/married-with-children
% https://www.irs.gov/newsroom/married-without-qualifying-children
if (bl_marital)
    % The benefit would start decreasing at a rate of $5 for every additional $100 in income
    fc_check_reduce = @(x) ((max(x - 150000,0))/100)*5;
end

% stimulus check function first check
fc_stimulus_check = @(x) max(0, fl_check_dollar - fc_check_reduce(x));

end

function ar_stimulus_check = snwi_stimulus_check_bush(...
    it_kids, bl_marital, fl_stimulus_adult, fl_stimulus_child, ar_income)

[mn_taxable_income, mn_tax_liability] = snw_tax_liability(ar_income);

kids = it_kids+1;
married = bl_marital+1;

ar_tax_liability = mn_tax_liability(:, married, it_kids+1);

% stimulus check function first check
% fl_stimulus_adult=600;
% fl_stimulus_child=300;
ar_stimulus_check = ...
    max(...
    max(300*(1+bl_marital), min(ar_tax_liability, fl_stimulus_adult*(1+bl_marital))) ...
    + it_kids*fl_stimulus_child ...
    - max(0, (ar_income' - 75000*(1+bl_marital)))*0.05, ...
    0);
end
