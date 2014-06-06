data exemplary;
	input trt mu;
	do obs=1 to 5;
		output;
	end;
datalines;
0 20
1 25
2 25
;

proc glimmix data=exemplary;
	class trt;
	model mu=trt / ddfm=KR;
	parms(9) / hold=1;
	contrast 'ctrl vs. exp' trt 2 -1 -1;
	contrast 'ctrl vs. exp 1' trt 1 -1 0;
	contrast 'ctrl vs. exp 2' trt 1 0 -1;
	contrast 'exp 1 vs. exp 2' trt 0 1 -1;
	lsmeans trt/diff cl;
	ods output tests3=F_overall contrasts=F_contrasts;
run;

data powerKR;
	set F_overall F_contrasts;
	nc_parm=numdf*Fvalue;
	alpha=0.05;
	F_Crit=Finv(1-alpha,numdf,dendf,0);
	Power=1-probF(F_crit,numdf,dendf,nc_parm);
run;

proc print data=powerKR;
run;
