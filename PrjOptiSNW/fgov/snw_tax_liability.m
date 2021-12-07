%% SNW_TAX_LIABILITY Given income find taxable income and tax liability
%    We can study the effects of the 2008 Tax Rebate. The Tax rebate is a
%    rebate based on how much tax was paid, so we need to know taxable
%    income and tax liability. These differ by income, household marital
%    status, and the count of children.
%
%    Given an array of pre-tax income values, we compute for from 0 to 4
%    kids and both married and unmarried taxable-income and tax-liability
%    at all points along the income array.
%
%    Deductions are from 2008 income
%    <https://www.irs.gov/pub/irs-prior/f1040--2008.pdf 2008 IRS f1040>.
%    Tax brackets from 2008 are here
%    <https://www.irs.gov/pub/irs-prior/i1040tt--2008.pdf 2008 IRS 1040tt>.
%
%    [MN_TAXABLE_INCOME, MN_TAX_LIABILITY] = SNW_TAX_LIABILITY(AR_INCOME,
%    BL_GRAPH) finds taxable income and tax-liability in 3d given income
%    array. 1st dimension is income, 2nd dimension is married or not (1st
%    column single, 2nd column married), and 3rd dimension is the number of
%    kids.
%
%    See also SNWX_TAX_LIABILITY
%

function [mn_taxable_income, mn_tax_liability] = snw_tax_liability(varargin)

if (~isempty(varargin))

    bl_graph = false;
    bl_verbose = false;
    if (length(varargin)==1)
        [ar_income] = varargin{:};
    elseif (length(varargin)==2)
        [ar_income, bl_graph] = varargin{:};
    elseif (length(varargin)==2)
        [ar_income, bl_graph, bl_verbose] = varargin{:};
    end

else

    clc;
    close all;
    clear all;

    bl_graph = true;
    bl_verbose = true;
    ar_income=linspace(0,400000,5)';
    ar_income = [20000, 30000, 50000, 70000];
    ar_income=linspace(20000,70000,100);
    ar_income=linspace(0,200000,100);
    ar_income=linspace(0,200000,5);
    %     ar_income = [20000];

end

%% A. Exemptions and Deductions
fl_exemption=3500; % exemption amount per household member

mt_deduction=zeros(2,5); % Marital-status and number of children-specific deduction
mt_deduction(1,1)=5450; % Single filers
mt_deduction(1,2:5)=8000; % Single filer with children
mt_deduction(2,:)=10900; % Married couples filing jointly

%% B. Taxable Income
% income = [70000];
mn_taxable_income = NaN(length(ar_income), 2, 5);

for y=1:length(ar_income)
    for m=1:2
        for k=0:4
            mn_taxable_income(y,m,k+1)=ar_income(y)-fl_exemption*m-fl_exemption*k-mt_deduction(m,k+1);
        end
    end
end

%% C. Tax Liability
mn_tax_liability=zeros(length(ar_income),2,5);

% Single filers with 0 children
for y=1:length(ar_income)
    if mn_taxable_income(y,1,1)<=8025
        mn_tax_liability(y,1,1)=0.1*mn_taxable_income(y,1,1);
    elseif ( mn_taxable_income(y,1,1)>8025 && mn_taxable_income(y,1,1)<=32550 )
        mn_tax_liability(y,1,1)=0.1*8025 + 0.15*( mn_taxable_income(y,1,1)-8025 );
    elseif ( mn_taxable_income(y,1,1)>32550 && mn_taxable_income(y,1,1)<=78850 )
        mn_tax_liability(y,1,1)=0.1*8025 + 0.15*(32550-8025) + 0.25*( mn_taxable_income(y,1,1)-32550 );
    elseif ( mn_taxable_income(y,1,1)>78850 && mn_taxable_income(y,1,1)<=164550 )
        mn_tax_liability(y,1,1)=0.1*8025 + 0.15*(32550-8025) + 0.25*(78850-32550) + 0.28*( mn_taxable_income(y,1,1)-78850 );
    elseif ( mn_taxable_income(y,1,1)>164550 && mn_taxable_income(y,1,1)<=357700 )
        mn_tax_liability(y,1,1)=0.1*8025 + 0.15*(32550-8025) + 0.25*(78850-32550) + 0.28*(164550-78850) + 0.33*( mn_taxable_income(y,1,1)-164550 );
    elseif ( mn_taxable_income(y,1,1)>357700 )
        mn_tax_liability(y,1,1)=0.1*8025 + 0.15*(32550-8025) + 0.25*(78850-32550) + 0.28*(164550-78850) + 0.33*(357700-164550) + 0.35*( mn_taxable_income(y,1,1)-357700 );
    end
end


% Single filers with 1-4 children
for k=2:5
    for y=1:length(ar_income)
        if mn_taxable_income(y,1,k)<=11450
            mn_tax_liability(y,1,k)=0.1*mn_taxable_income(y,1,k);
        elseif ( mn_taxable_income(y,1,k)>11450 && mn_taxable_income(y,1,k)<=43650 )
            mn_tax_liability(y,1,k)=0.1*11450 + 0.15*( mn_taxable_income(y,1,k)-11450 );
        elseif ( mn_taxable_income(y,1,k)>43650 && mn_taxable_income(y,1,k)<=112650 )
            mn_tax_liability(y,1,k)=0.1*11450 + 0.15*(43650-11450) + 0.25*( mn_taxable_income(y,1,k)-43650 );
        elseif ( mn_taxable_income(y,1,k)>112650 && mn_taxable_income(y,1,k)<=182400 )
            mn_tax_liability(y,1,k)=0.1*11450 + 0.15*(43650-11450) + 0.25*(112650-43650) + 0.28*( mn_taxable_income(y,1,k)-112650 );
        elseif ( mn_taxable_income(y,1,k)>182400 && mn_taxable_income(y,1,k)<=357700 )
            mn_tax_liability(y,1,k)=0.1*11450 + 0.15*(43650-11450) + 0.25*(112650-43650) + 0.28*(182400-112650) + 0.33*( mn_taxable_income(y,1,k)-182400 );
        elseif ( mn_taxable_income(y,1,k)>357700 )
            mn_tax_liability(y,1,k)=0.1*11450 + 0.15*(43650-11450) + 0.25*(112650-43650) + 0.28*(182400-112650) + 0.33*(357700-182400) + 0.35*( mn_taxable_income(y,1,k)-357700 );
        end
    end
end


% Married filers
for k=1:5
    for y=1:length(ar_income)
        if mn_taxable_income(y,2,k)<=16050
            mn_tax_liability(y,2,k)=0.1*mn_taxable_income(y,2,k);
        elseif ( mn_taxable_income(y,2,k)>16050 && mn_taxable_income(y,2,k)<=65100 )
            mn_tax_liability(y,2,k)=0.1*16050 + 0.15*( mn_taxable_income(y,2,k)-16050 );
        elseif ( mn_taxable_income(y,2,k)>65100 && mn_taxable_income(y,2,k)<=131450 )
            mn_tax_liability(y,2,k)=0.1*16050 + 0.15*(65100-16050) + 0.25*( mn_taxable_income(y,2,k)-65100 );
        elseif ( mn_taxable_income(y,2,k)>131450 && mn_taxable_income(y,2,k)<=200300 )
            mn_tax_liability(y,2,k)=0.1*16050 + 0.15*(65100-16050) + 0.25*(131450-65100) + 0.28*( mn_taxable_income(y,2,k)-131450 );
        elseif ( mn_taxable_income(y,2,k)>200300 && mn_taxable_income(y,2,k)<=357700 )
            mn_tax_liability(y,2,k)=0.1*16050 + 0.15*(65100-16050) + 0.25*(131450-65100) + 0.28*(200300-131450) + 0.33*( mn_taxable_income(y,2,k)-200300 );
        elseif ( mn_taxable_income(y,2,k)>357700 )
            mn_tax_liability(y,2,k)=0.1*16050 + 0.15*(65100-16050) + 0.25*(131450-65100) + 0.28*(200300-131450) + 0.33*(357700-200300) + 0.35*( mn_taxable_income(y,2,k)-357700 );
        end
    end
end

for y=1:length(ar_income)
    for m=1:2
        for k=0:4
            mn_tax_liability(y,m,k+1)=max(mn_tax_liability(y,m,k+1) , 0);
        end
    end
end

%% D. Show Results and Visualization

if (bl_verbose)
    disp('mn_taxable_income');
    disp(mn_taxable_income);
    disp('mn_tax_liability');
    disp(mn_tax_liability);
end

% stimulus = max(max(300*(1+married), min(tax_liablity, 600*(1+married))) +  kids*300 - max(0, (ar_income - 75000*(1+married))*0.05), 0)
if (bl_graph)
    % Figure 1
    % maybe two graphs: one with taxable income on the vertical axis and income on the horizontal axis (for the various household types)

    mt_value = reshape(mn_taxable_income(:,1,:), [size(mn_taxable_income,1), size(mn_taxable_income,3)]);
    mt_value = mt_value';
    ar_row_grid = ["0 kids" "1 kid", "2 kids", "3 kids", "4 kids"];
    ar_col_grid = ar_income;
    mp_support_graph = containers.Map('KeyType', 'char', 'ValueType', 'any');
    mp_support_graph('cl_st_graph_title') = {'Income and Taxable Income (unmarried)'};
    mp_support_graph('cl_st_ytitle') = {'Taxable Income'};
    mp_support_graph('cl_st_xtitle') = {'Income'};
    mp_support_graph('bl_graph_logy') = false; % do not log
    ff_graph_grid(mt_value, ar_row_grid, ar_col_grid, mp_support_graph);

    mt_value = reshape(mn_taxable_income(:,2,:), [size(mn_taxable_income,1), size(mn_taxable_income,3)]);
    mt_value = mt_value';
    ar_row_grid = ["0 kids" "1 kid", "2 kids", "3 kids", "4 kids"];
    ar_col_grid = ar_income;
    mp_support_graph = containers.Map('KeyType', 'char', 'ValueType', 'any');
    mp_support_graph('cl_st_graph_title') = {'Income and Taxable Income (married)'};
    mp_support_graph('cl_st_ytitle') = {'Taxable Income'};
    mp_support_graph('cl_st_xtitle') = {'Income'};
    mp_support_graph('bl_graph_logy') = false; % do not log
    ff_graph_grid(mt_value, ar_row_grid, ar_col_grid, mp_support_graph);

    % Figure 2
    % and one graph with tax liability on the vertical axis and taxable income on the horizontal axis (for the various household types)

    mt_value = reshape(mn_tax_liability(:,1,:), [size(mn_tax_liability,1), size(mn_tax_liability,3)]);
    mt_value = mt_value';
    ar_row_grid = ["0 kids" "1 kid", "2 kids", "3 kids", "4 kids"];
    ar_col_grid = ar_income;
    mp_support_graph = containers.Map('KeyType', 'char', 'ValueType', 'any');
    mp_support_graph('cl_st_graph_title') = {'Income and Taxable Liability (unmarried)'};
    mp_support_graph('cl_st_ytitle') = {'Taxable Liability'};
    mp_support_graph('cl_st_xtitle') = {'Income'};
    mp_support_graph('bl_graph_logy') = false; % do not log
    ff_graph_grid(mt_value, ar_row_grid, ar_col_grid, mp_support_graph);

    mt_value = reshape(mn_tax_liability(:,2,:), [size(mn_tax_liability,1), size(mn_tax_liability,3)]);
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
