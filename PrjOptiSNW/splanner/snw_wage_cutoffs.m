function cutoffs=snw_wage_cutoffs(Phi_true, ...
    theta, epsilon, eta_H_grid, n_agrid, n_etagrid, n_educgrid, n_marriedgrid, n_kidsgrid, jret)

% global theta epsilon eta_H_grid n_agrid n_etagrid n_educgrid n_marriedgrid n_kidsgrid jret

Phi=zeros((jret-1)*n_agrid*n_etagrid*n_educgrid*n_marriedgrid*n_kidsgrid,2);
counter=0;

for j=1:(jret-1) % Age
   for a=1:n_agrid % Assets
       for eta=1:n_etagrid % Productivity
           for educ=1:n_educgrid % Educational level
               for married=1:n_marriedgrid % Marital status
                   for kids=1:n_kidsgrid % No. of kids
                       
                       counter=counter+1;
                                              
                       Phi(counter,1)=epsilon(j,educ)*theta*exp(eta_H_grid(eta));
                       Phi(counter,2)=Phi_true(j,a,eta,educ,married,kids);
                       
                   end
               end
           end
       end
   end
end

% Sort matrix
[~,idx]=sort(Phi(:,1));
sortedPhi=Phi(idx,:);

% Normalize to sum to 1
sortedPhi_sum=sum(sortedPhi(:,2));
sortedPhi(:,2)=sortedPhi(:,2)/sortedPhi_sum;
% Cumulative mass
sortedPhi(1,3)=sortedPhi(1,2);
for i=2:length(sortedPhi)
    sortedPhi(i,3)=sortedPhi(i-1,3)+sortedPhi(i,2);
end

clear sortedPhi_sum

ind_aux(1)=find(sortedPhi(:,3)<=0.2,1,'last');
ind_aux(2)=find(sortedPhi(:,3)<=0.4,1,'last');
ind_aux(3)=find(sortedPhi(:,3)<=0.6,1,'last');
ind_aux(4)=find(sortedPhi(:,3)<=0.8,1,'last');

cutoffs(1)=sortedPhi(ind_aux(1),1);
cutoffs(2)=sortedPhi(ind_aux(2),1);
cutoffs(3)=sortedPhi(ind_aux(3),1);
cutoffs(4)=sortedPhi(ind_aux(4),1);

name='Wage quintile cutoffs=';
name2=[name,num2str(cutoffs(:)')];
disp(name2);

end