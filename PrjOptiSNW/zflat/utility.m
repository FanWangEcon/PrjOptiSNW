function F=utility(cons,m,k)

global gamma cons_allocation_rule

hh_size=m+k-1; % m=1 if single; m=2 if married; k=1 if 0 children

if cons_allocation_rule==1 % Uniform consumption allocation rule
    
    if gamma==1
        if cons<=0
            F=-10E30;
        else
            F=log(cons/hh_size);
        end
    else
        if cons<=0
            F=-10E30;
        else
            F=(((cons/hh_size).^(1-gamma))-1)/(1-gamma);
        end
    end
    
elseif cons_allocation_rule==2 % Square root consumption allocation rule
    
    if gamma==1
        if cons<=0
            F=-10E30;
        else
            F=log(cons/sqrt(hh_size));
        end
    else
        if cons<=0
            F=-10E30;
        else
            F=(((cons/sqrt(hh_size)).^(1-gamma))-1)/(1-gamma);
        end
    end
    
end

end