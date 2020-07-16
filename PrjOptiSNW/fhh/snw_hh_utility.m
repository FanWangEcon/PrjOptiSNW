function F=snw_hh_utility(cons,m,k,gamma, cons_allocation_rule)

% global gamma cons_allocation_rule

hh_size=m+k-1; % m=1 if single; m=2 if married; k=1 if 0 children

if cons_allocation_rule==0 % Household size does not matter

    hh_power=0;
    
%     if gamma==1
%         if cons<=0
%             F=-10E30;
%         else
%             F=log(cons/hh_size);
%         end
%     else
%         if cons<=0
%             F=-10E30;
%         else
%             F=(((cons/hh_size).^(1-gamma))-1)/(1-gamma);
%         end
%     end
    
else
 
    hh_power=1/cons_allocation_rule;    
    
%     if gamma==1
%         if cons<=0
%             F=-10E30;
%         else
%             F=log(cons/sqrt(hh_size));
%         end
%     else
%         if cons<=0
%             F=-10E30;
%         else
%             F=(((cons/sqrt(hh_size)).^(1-gamma))-1)/(1-gamma);
%         end
%     end
    
% elseif cons_allocation_rule==0 % Square root consumption allocation rule
%     
%     hh_power=0;
    
end

if gamma==1
    if cons<=0
        F=-10E30;
    else
        F=log(cons/hh_size^hh_power);
    end
else
    if cons<=0
        F=-10E30;
    else
        F=(((cons/hh_size^hh_power).^(1-gamma))-1)/(1-gamma);
    end
end

end