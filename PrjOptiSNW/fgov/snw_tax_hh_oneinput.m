%% SNW_TAX_HH Tax levels given household head and spousal incomes
%    Level of taxes
%

function F=snw_tax_hh_oneinput(varargin)

if (~isempty(varargin))
    
    a0=0.258;
    a1=0.768;
    bl_snw_tax_hh_verbose = false;
    
    if (length(varargin)==2)
        [inc_plus_spouse_inc, a2] = varargin{:};
    elseif (length(varargin)==3)
        [inc_plus_spouse_inc, a2, a0] = varargin{:};
    elseif (length(varargin)==4)
        [inc_plus_spouse_inc, a2, a0, a1] = varargin{:};
    elseif (length(varargin)==5)
        [inc_plus_spouse_inc, a2, a0, a1, bl_snw_tax_hh_verbose] = varargin{:};
    end
    
    clear varargin
else
    
    clc;
    inc = 1;
    spouse_inc = 0;
    a0=0.258;
    a1=0.768;    
    a2 = 1.5286;
    bl_snw_tax_hh_verbose = true;
    
end

F=a0*(inc_plus_spouse_inc-(((inc_plus_spouse_inc.^(-a1))+a2).^(-1/a1)));

if(bl_snw_tax_hh_verbose)
    
    figure();
    close all;
    hold on;
    x_low = 0;
    x_high = 5;
    a2 = 1.5286;
    fplot(@(x) a0*(x-(((x.^(-a1))+a2).^(-1/a1)))/x, [x_low, x_high]);
    a2 = 12.7176;
    fplot(@(x) a0*(x-(((x.^(-a1))+a2).^(-1/a1)))/x, [x_low, x_high]);
    xlabel('income');
    ylabel('tax rate');
    title(['Tax Rate and Income'])
    legend('a2=1.5286','a2=12.7176');
    grid on
    
    figure();
    hold on;
    a2 = 1.5286;
    a0 = 0.258;
    fplot(@(x) a0*(x-(((x.^(-a1))+a2).^(-1/a1)))/x, [x_low, x_high]);
    a0 = 0.31739;
    fplot(@(x) a0*(x-(((x.^(-a1))+a2).^(-1/a1)))/x, [x_low, x_high]);
    xlabel('income');
    ylabel('tax rate');
    title(['Tax Rate and Income'])
    legend('a2 = 1.5286, a0 = 0.258','a2 = 1.5286, a0 = 0.31739');
    grid on    
    
end

end