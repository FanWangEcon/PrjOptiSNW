function [factor, cutoffs]=pi_unemp_calibration(Phi_true,pi_j,pi_w)

global theta epsilon eta_grid n_agrid n_etagrid n_educgrid n_marriedgrid n_kidsgrid

Phi=zeros(24*n_agrid*n_etagrid*n_educgrid*n_marriedgrid*n_kidsgrid,2);
counter=0;

for j=1:24 % Age
   for a=1:n_agrid % Assets
       for eta=1:n_etagrid % Productivity
           for educ=1:n_educgrid % Educational level
               for married=1:n_marriedgrid % Marital status
                   for kids=1:n_kidsgrid % No. of kids
                       
                       counter=counter+1;
                                              
                       Phi(counter,1)=epsilon(j,educ)*theta*exp(eta_grid(eta));
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

clear sortedPhi_sum ind_aux name name2

% Find distribution over wage quintiles given age group
Phi_aux=zeros(24,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid,5);

for j=1:24 % Age
   for a=1:n_agrid % Assets
       for eta=1:n_etagrid % Productivity
           for educ=1:n_educgrid % Educational level
               for married=1:n_marriedgrid % Marital status
                   for kids=1:n_kidsgrid % No. of kids
                       
                       wages=epsilon(j,educ)*theta*exp(eta_grid(eta));
                       if wages<=cutoffs(1)
                           Phi_aux(j,a,eta,educ,married,kids,1)=Phi_aux(j,a,eta,educ,married,kids,1)+Phi_true(j,a,eta,educ,married,kids);
                       elseif wages>cutoffs(1) && wages<=cutoffs(2)
                           Phi_aux(j,a,eta,educ,married,kids,2)=Phi_aux(j,a,eta,educ,married,kids,2)+Phi_true(j,a,eta,educ,married,kids);
                       elseif wages>cutoffs(2) && wages<=cutoffs(3)
                           Phi_aux(j,a,eta,educ,married,kids,3)=Phi_aux(j,a,eta,educ,married,kids,3)+Phi_true(j,a,eta,educ,married,kids);
                       elseif wages>cutoffs(3) && wages<=cutoffs(4)
                           Phi_aux(j,a,eta,educ,married,kids,4)=Phi_aux(j,a,eta,educ,married,kids,4)+Phi_true(j,a,eta,educ,married,kids);
                       elseif wages>cutoffs(4)
                           Phi_aux(j,a,eta,educ,married,kids,5)=Phi_aux(j,a,eta,educ,married,kids,5)+Phi_true(j,a,eta,educ,married,kids);
                       end
                                              
                   end
               end
           end
       end
   end
end

% Normalize mass to sum to 1 for each age group
Phi_aux_j1=sum(sum(sum(sum(sum(sum(sum(Phi_aux(1:7,:,:,:,:,:,:))))))));
Phi_aux_j2=sum(sum(sum(sum(sum(sum(sum(Phi_aux(8:12,:,:,:,:,:,:))))))));
Phi_aux_j3=sum(sum(sum(sum(sum(sum(sum(Phi_aux(13:17,:,:,:,:,:,:))))))));
Phi_aux_j4=sum(sum(sum(sum(sum(sum(sum(Phi_aux(18:22,:,:,:,:,:,:))))))));
Phi_aux_j5=sum(sum(sum(sum(sum(sum(sum(Phi_aux(23:24,:,:,:,:,:,:))))))));

Phi_aux_j_norm=zeros(24,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid,5);
Phi_aux_j_norm(1:7,:,:,:,:,:,:)=Phi_aux(1:7,:,:,:,:,:,:)/Phi_aux_j1;
Phi_aux_j_norm(8:12,:,:,:,:,:,:)=Phi_aux(8:12,:,:,:,:,:,:)/Phi_aux_j2;
Phi_aux_j_norm(13:17,:,:,:,:,:,:)=Phi_aux(13:17,:,:,:,:,:,:)/Phi_aux_j3;
Phi_aux_j_norm(18:22,:,:,:,:,:,:)=Phi_aux(18:22,:,:,:,:,:,:)/Phi_aux_j4;
Phi_aux_j_norm(23:24,:,:,:,:,:,:)=Phi_aux(23:24,:,:,:,:,:,:)/Phi_aux_j5;

pi_w_1=NaN(1,5);
pi_w_2=NaN(1,5);
pi_w_3=NaN(1,5);
pi_w_4=NaN(1,5);
pi_w_5=NaN(1,5);

for i=1:5
    pi_w_1(i)=sum(sum(sum(sum(sum(sum(Phi_aux_j_norm(1:7,:,:,:,:,:,i)))))));
    pi_w_2(i)=sum(sum(sum(sum(sum(sum(Phi_aux_j_norm(8:12,:,:,:,:,:,i)))))));
    pi_w_3(i)=sum(sum(sum(sum(sum(sum(Phi_aux_j_norm(13:17,:,:,:,:,:,i)))))));
    pi_w_4(i)=sum(sum(sum(sum(sum(sum(Phi_aux_j_norm(18:22,:,:,:,:,:,i)))))));
    pi_w_5(i)=sum(sum(sum(sum(sum(sum(Phi_aux_j_norm(23:24,:,:,:,:,:,i)))))));
end

check1=0;
check2=0;
check3=0;
check4=0;
check5=0;

for i=1:5
   check1=check1+(pi_w(i)*pi_w_1(i));
   check2=check2+(pi_w(i)*pi_w_2(i));
   check3=check3+(pi_w(i)*pi_w_3(i));
   check4=check4+(pi_w(i)*pi_w_4(i));
   check5=check5+(pi_w(i)*pi_w_5(i));
end

factor(1)=pi_j(1)/check1;
factor(2)=pi_j(2)/check2;
factor(3)=pi_j(3)/check3;
factor(4)=pi_j(4)/check4;
factor(5)=pi_j(5)/check5;

end