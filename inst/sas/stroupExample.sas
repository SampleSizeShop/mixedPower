data exemplary;
	input id trt time mu;
datalines;
1 0 1 20
1 0 2 20
1 0 3 20
2 0 1 20
2 0 2 20
2 0 3 20
3 0 1 20
3 0 2 20
3 0 3 20
4 0 1 20
4 0 2 20
4 0 3 20
5 0 1 20
5 0 2 20
5 0 3 20
6 1 1 20
6 1 2 23
6 1 3 25
7 1 1 20
7 1 2 23
7 1 3 25
8 1 1 20
8 1 2 23
8 1 3 25
9 1 1 20
9 1 2 23
9 1 3 25
10 1 1 20
10 1 2 23 
10 1 3 25
11 2 1 20
11 2 2 23 
11 2 3 25
12 2 1 20
12 2 2 23 
12 2 3 25
13 2 1 20
13 2 2 23 
13 2 3 25
14 2 1 20
14 2 2 23 
14 2 3 25
15 2 1 20
15 2 2 23 
15 2 3 25
;
run;

data exemplary;
	set exemplary;
	trt0 = (trt = 0);
	trt1 = (trt = 1);
	trt2 = (trt = 2);
	time1 = (time = 1);
	time2 = (time = 2);
	time3 = (time = 3);
	trt0_time1 = trt0 * time1;
	trt0_time2 = trt0 * time2;
	trt0_time3 = trt0 * time3;
	trt1_time1 = trt1 * time1;
	trt1_time2 = trt1 * time2;
	trt1_time3 = trt1 * time3;
	trt2_time1 = trt2 * time1;
	trt2_time2 = trt2 * time2;
	trt2_time3 = trt2 * time3;
run;

proc glimmix data=exemplary;
	class id;
	model mu = 
		trt0_time1 trt0_time2 trt0_time3
		trt1_time1 trt1_time2 trt1_time3
		trt2_time1 trt2_time2 trt2_time3

		/ noint ddfm=KR;
	random _residual_ / subject=id type=CS;
	parms (2) (0.7) / hold=1,2 noiter;

	contrast 'time by trt' 
		trt0_time1 1 trt0_time2 -1 trt0_time3 0
		trt1_time1 -1 trt1_time2 1 trt1_time3 0
		trt2_time1 0 trt2_time2 0 trt2_time3 0 ,

		trt0_time1 1 trt0_time2 0 trt0_time3 -1
		trt1_time1 -1 trt1_time2 0 trt1_time3 1
		trt2_time1 0 trt2_time2 0 trt2_time3 0,

		trt0_time1 1 trt0_time2 -1 trt0_time3 0
		trt1_time1 0 trt1_time2 0 trt1_time3 0
		trt2_time1 -1 trt2_time2 1 trt2_time3 0,

		trt0_time1 1 trt0_time2 0 trt0_time3 -1
		trt1_time1 0 trt1_time2 0 trt1_time3 0
		trt2_time1 -1 trt2_time2 0 trt2_time3 1
		;
	ods output tests3=F_overall contrasts=F_contrasts;
run;

/*
proc glimmix data=exemplary;
	class trt;
	model mu=trt / ddfm=BW;
	parms(9) / hold=1;
	contrast 'ctrl vs. exp' trt 2 -1 -1;
	contrast 'ctrl vs. exp 1' trt 1 -1 0;
	contrast 'ctrl vs. exp 2' trt 1 0 -1;
	contrast 'exp 1 vs. exp 2' trt 0 1 -1;
	lsmeans trt/diff cl;
	ods output tests3=F_overall contrasts=F_contrasts;
run;
*/
data powerKR;
	set F_contrasts;
	nc_parm=numdf*Fvalue;
	alpha=0.05;
	F_Crit=Finv(1-alpha,numdf,dendf,0);
	Power=1-probF(F_crit,numdf,dendf,nc_parm);
run;

proc print data=powerKR;
run;
