%% SNW_TAX_HH Tax levels given household head and spousal incomes
%    Level of taxes
%

function F=snw_tax_hh(varargin)

if (~isempty(varargin))
    
    a0=0.258;
    a1=0.768;
    bl_snw_tax_hh_verbose = false;
    
    if (length(varargin)==3)
        [inc, spouse_inc, a2] = varargin{:};
    elseif (length(varargin)==4)
        [inc, spouse_inc, a2, a0] = varargin{:};
    elseif (length(varargin)==5)
        [inc, spouse_inc, a2, a0, a1] = varargin{:};
    elseif (length(varargin)==6)
        [inc, spouse_inc, a2, a0, a1, bl_snw_tax_hh_verbose] = varargin{:};
    end
    
%     clear varargin
else
    
    clc;
    inc = 1;
    spouse_inc = 0;
    a0=0.258;
    a1=0.768;    
    a2 = 1.5286;
    bl_snw_tax_hh_verbose = true;
    
end

inc_fam=inc+spouse_inc;
% clear inc spouse_inc
F=a0*(inc_fam-(((inc_fam.^(-a1))+a2).^(-1/a1)));

if(bl_snw_tax_hh_verbose)
    
    figure();
    close all;
    hold on;
    x_low = 0;
    x_high = 5;
    a2 = 1.5286;
    fplot(@(x) a0*(x-(((x.^(-a1))+a2).^(-1/a1)))/x, [x_low, x_high]);
    a2 = 1.6213;
    fplot(@(x) a0*(x-(((x.^(-a1))+a2).^(-1/a1)))/x, [x_low, x_high]);
    xlabel('income');
    ylabel('tax rate');
    title(['Tax Rate and Income'])
    legend('prior-0.97beta-a2=1.5286','in-progress-new-a2=1.6213');
    grid on
    
%     figure();
%     hold on;
%     a2 = 1.5286;
%     a0 = 0.258;
%     fplot(@(x) a0*(x-(((x.^(-a1))+a2).^(-1/a1)))/x, [x_low, x_high]);
%     a0 = 0.31739;
%     fplot(@(x) a0*(x-(((x.^(-a1))+a2).^(-1/a1)))/x, [x_low, x_high]);
%     xlabel('income');
%     ylabel('tax rate');
%     title(['Tax Rate and Income'])
%     legend('a2 = 1.5286, a0 = 0.258','a2 = 1.5286, a0 = 0.31739');
%     grid on    
    
end

end