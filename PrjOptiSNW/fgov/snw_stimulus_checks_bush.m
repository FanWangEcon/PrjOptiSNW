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

exemption=3500; % exemption amount per household member

deduction=zeros(2,5); % Marital-status and number of children-specific deduction
deduction(1,1)=5450; % Single filers
deduction(1,2:5)=8000; % Single filer with children
deduction(2,:)=10900; % Married couples filing jointly

% ar_income=linspace(0,400000,5)';
% ar_income = [20000, 30000, 50000, 70000];
% ar_income=linspace(20000,70000,100);
% ar_income=linspace(0,200000,100);
% income = [20000];

% income = [70000];

taxable_income = NaN(length(ar_income), 2, 5);
for y=1:length(ar_income)
    for m=1:2
        for k=0:4
            taxable_income(y,m,k+1)=ar_income(y)-exemption*m-exemption*k-deduction(m,k+1);
        end
    end
end

tax_liability=zeros(length(ar_income),2,5);

% Single filers with 0 children
for y=1:length(ar_income)
    if taxable_income(y,1,1)<=8025
        tax_liability(y,1,1)=0.1*taxable_income(y,1,1);
    elseif ( taxable_income(y,1,1)>8025 && taxable_income(y,1,1)<=32550 )
        tax_liability(y,1,1)=0.1*8025 + 0.15*( taxable_income(y,1,1)-8025 );
    elseif ( taxable_income(y,1,1)>32550 && taxable_income(y,1,1)<=78850 )
        tax_liability(y,1,1)=0.1*8025 + 0.15*(32550-8025) + 0.25*( taxable_income(y,1,1)-32550 );
    elseif ( taxable_income(y,1,1)>78850 && taxable_income(y,1,1)<=164550 )
        tax_liability(y,1,1)=0.1*8025 + 0.15*(32550-8025) + 0.25*(78850-32550) + 0.28*( taxable_income(y,1,1)-78850 );
    elseif ( taxable_income(y,1,1)>164550 && taxable_income(y,1,1)<=357700 )
        tax_liability(y,1,1)=0.1*8025 + 0.15*(32550-8025) + 0.25*(78850-32550) + 0.28*(164550-78850) + 0.33*( taxable_income(y,1,1)-164550 );
    elseif ( taxable_income(y,1,1)>357700 )
        tax_liability(y,1,1)=0.1*8025 + 0.15*(32550-8025) + 0.25*(78850-32550) + 0.28*(164550-78850) + 0.33*(357700-164550) + 0.35*( taxable_income(y,1,1)-357700 );
    end
end


% Single filers with 1-4 children
for k=2:5
    for y=1:length(ar_income)
        if taxable_income(y,1,k)<=11450
            tax_liability(y,1,k)=0.1*taxable_income(y,1,k);
        elseif ( taxable_income(y,1,k)>11450 && taxable_income(y,1,k)<=43650 )
            tax_liability(y,1,k)=0.1*11450 + 0.15*( taxable_income(y,1,k)-11450 );
        elseif ( taxable_income(y,1,k)>43650 && taxable_income(y,1,k)<=112650 )
            tax_liability(y,1,k)=0.1*11450 + 0.15*(43650-11450) + 0.25*( taxable_income(y,1,k)-43650 );
        elseif ( taxable_income(y,1,k)>112650 && taxable_income(y,1,k)<=182400 )
            tax_liability(y,1,k)=0.1*11450 + 0.15*(43650-11450) + 0.25*(112650-43650) + 0.28*( taxable_income(y,1,k)-112650 );
        elseif ( taxable_income(y,1,k)>182400 && taxable_income(y,1,k)<=357700 )
            tax_liability(y,1,k)=0.1*11450 + 0.15*(43650-11450) + 0.25*(112650-43650) + 0.28*(182400-112650) + 0.33*( taxable_income(y,1,k)-182400 );
        elseif ( taxable_income(y,1,k)>357700 )
            tax_liability(y,1,k)=0.1*11450 + 0.15*(43650-11450) + 0.25*(112650-43650) + 0.28*(182400-112650) + 0.33*(357700-182400) + 0.35*( taxable_income(y,1,k)-357700 );
        end
    end
end


% Married filers
for k=1:5
    for y=1:length(ar_income)
        if taxable_income(y,2,k)<=16050
            tax_liability(y,2,k)=0.1*taxable_income(y,2,k);
        elseif ( taxable_income(y,2,k)>16050 && taxable_income(y,2,k)<=65100 )
            tax_liability(y,2,k)=0.1*16050 + 0.15*( taxable_income(y,2,k)-16050 );
        elseif ( taxable_income(y,2,k)>65100 && taxable_income(y,2,k)<=131450 )
            tax_liability(y,2,k)=0.1*16050 + 0.15*(65100-16050) + 0.25*( taxable_income(y,2,k)-65100 );
        elseif ( taxable_income(y,2,k)>131450 && taxable_income(y,2,k)<=200300 )
            tax_liability(y,2,k)=0.1*16050 + 0.15*(65100-16050) + 0.25*(131450-65100) + 0.28*( taxable_income(y,2,k)-131450 );
        elseif ( taxable_income(y,2,k)>200300 && taxable_income(y,2,k)<=357700 )
            tax_liability(y,2,k)=0.1*16050 + 0.15*(65100-16050) + 0.25*(131450-65100) + 0.28*(200300-131450) + 0.33*( taxable_income(y,2,k)-200300 );
        elseif ( taxable_income(y,2,k)>357700 )
            tax_liability(y,2,k)=0.1*16050 + 0.15*(65100-16050) + 0.25*(131450-65100) + 0.28*(200300-131450) + 0.33*(357700-200300) + 0.35*( taxable_income(y,2,k)-357700 );
        end
    end
end

for y=1:length(ar_income)
    for m=1:2
        for k=0:4
            tax_liability(y,m,k+1)=max(tax_liability(y,m,k+1) , 0);
        end
    end
end

% disp(taxable_income);
% 
% disp(tax_liability);


kids = it_kids+1;
married = bl_marital+1;

ar_tax_liability = tax_liability(:, married, it_kids+1);

% stimulus check function first check
% fl_stimulus_adult=600;
% fl_stimulus_child=300;
ar_stimulus_check = ...
    max(...
    max(300*(1+bl_marital), min(ar_tax_liability, fl_stimulus_adult*(1+bl_marital))) ...
    + it_kids*fl_stimulus_child ...
    - max(0, (ar_income' - 75000*(1+bl_marital)))*0.05, ...
    0);

%% Figure 1
% maybe two graphs: one with taxable income on the vertical axis and income on the horizontal axis (for the various household types)

bl_graph = false;
if (bl_graph)
    mt_value = reshape(taxable_income(:,1,:), [size(taxable_income,1), size(taxable_income,3)]);
    mt_value = mt_value';
    ar_row_grid = ["0 kids" "1 kid", "2 kids", "3 kids", "4 kids"];
    ar_col_grid = ar_income;
    mp_support_graph = containers.Map('KeyType', 'char', 'ValueType', 'any');
    mp_support_graph('cl_st_graph_title') = {'Income and Taxable Income (unmarried)'};
    mp_support_graph('cl_st_ytitle') = {'Taxable Income'};
    mp_support_graph('cl_st_xtitle') = {'Income'};
    mp_support_graph('bl_graph_logy') = false; % do not log
    ff_graph_grid(mt_value, ar_row_grid, ar_col_grid, mp_support_graph);
    
    mt_value = reshape(taxable_income(:,2,:), [size(taxable_income,1), size(taxable_income,3)]);
    mt_value = mt_value';
    ar_row_grid = ["0 kids" "1 kid", "2 kids", "3 kids", "4 kids"];
    ar_col_grid = ar_income;
    mp_support_graph = containers.Map('KeyType', 'char', 'ValueType', 'any');
    mp_support_graph('cl_st_graph_title') = {'Income and Taxable Income (married)'};
    mp_support_graph('cl_st_ytitle') = {'Taxable Income'};
    mp_support_graph('cl_st_xtitle') = {'Income'};
    mp_support_graph('bl_graph_logy') = false; % do not log
    ff_graph_grid(mt_value, ar_row_grid, ar_col_grid, mp_support_graph);
    
    %% Figure 2
    % and one graph with tax liability on the vertical axis and taxable income on the horizontal axis (for the various household types)
    
    mt_value = reshape(tax_liability(:,1,:), [size(tax_liability,1), size(tax_liability,3)]);
    mt_value = mt_value';
    ar_row_grid = ["0 kids" "1 kid", "2 kids", "3 kids", "4 kids"];
    ar_col_grid = ar_income;
    mp_support_graph = containers.Map('KeyType', 'char', 'ValueType', 'any');
    mp_support_graph('cl_st_graph_title') = {'Income and Taxable Liability (unmarried)'};
    mp_support_graph('cl_st_ytitle') = {'Taxable Liability'};
    mp_support_graph('cl_st_xtitle') = {'Income'};
    mp_support_graph('bl_graph_logy') = false; % do not log
    ff_graph_grid(mt_value, ar_row_grid, ar_col_grid, mp_support_graph);
    
    mt_value = reshape(tax_liability(:,2,:), [size(tax_liability,1), size(tax_liability,3)]);
    mt_value = mt_value';
    ar_row_grid = ["0 kids" "1 kid", "2 kids", "3 kids", "4 kids"];
    ar_col_grid = ar_income;
    mp_support_graph = containers.Map('KeyType', 'char', 'ValueType', 'any');
    mp_support_graph('cl_st_graph_title') = {'Income and Taxable Liability (married)'};
    mp_support_graph('cl_st_ytitle') = {'Taxable Liability'};
    mp_support_graph('cl_st_xtitle') = {'Income'};
    mp_support_graph('bl_graph_logy') = false; % do not log
    ff_graph_grid(mt_value, ar_row_grid, ar_col_grid, mp_support_graph);
end
end