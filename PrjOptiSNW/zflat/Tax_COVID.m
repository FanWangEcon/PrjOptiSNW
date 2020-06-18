function F=Tax_COVID(inc,spouse_inc,a2_COVID)

a0=0.258;
a1=0.768;

inc_fam=inc+spouse_inc;

F=a0*(inc_fam-(((inc_fam^(-a1))+a2_COVID)^(-1/a1)));

end