function F=Tax(inc,spouse_inc)

global a2

a0=0.258;
a1=0.768;

inc_fam=inc+spouse_inc;

F=a0*(inc_fam-(((inc_fam^(-a1))+a2)^(-1/a1)));

end