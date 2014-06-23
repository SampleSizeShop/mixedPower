  * create exemplary data for complete cases in a 2 group, 5 repeated measures design;
  data exemplaryData;
	do grp = 1 to 2;
	  do isu = 1 to 50; 
	    do rep = 1 to 5;
			y = 2.5 * (grp = 1 & rep = 1);
			subjectID = ((grp-1) * 50) + isu;
			trt1_rep1 = (grp = 1 & rep = 1);
			trt1_rep2 = (grp = 1 & rep = 2);
			trt1_rep3 = (grp = 1 & rep = 3);
			trt1_rep4 = (grp = 1 & rep = 4);
			trt1_rep5 = (grp = 1 & rep = 5);
  			trt2_rep1 = (grp = 2 & rep = 1); 
			trt2_rep2 = (grp = 2 & rep = 2);
			trt2_rep3 = (grp = 2 & rep = 3);
			trt2_rep4 = (grp = 2 & rep = 4);
			trt2_rep5 = (grp = 2 & rep = 5);
			output;
		end;
	  end;
    end;
	drop isu grp rep;
  run;

proc glimmix data=exemplaryData;
	class subjectID;
    	model y = trt1_rep1 trt1_rep2 trt1_rep3 trt1_rep4 trt1_rep5
  				trt2_rep1 trt2_rep2 trt2_rep3 trt2_rep4 trt2_rep5 / noint solution ddfm=kr;
  		random _residual_ / subject=subjectID type=CS;
      parms (1) (0.4) / hold=1,2;
  		contrast "time by treatment"
  			trt1_rep1 1 trt1_rep2 -1 trt2_rep1 -1 trt2_rep2 1,
  			trt1_rep1 1 trt1_rep3 -1 trt2_rep1 -1 trt2_rep3 1,
  			trt1_rep1 1 trt1_rep4 -1 trt2_rep1 -1 trt2_rep4 1,
  			trt1_rep1 1 trt1_rep5 -1 trt2_rep1 -1 trt2_rep5 1;
	  ods output contrasts=tmpContrasts;
run;

  data power;
    set tmpContrasts;
    Noncen = NumDF * FValue;
    Alpha = 0.05;
    FCrit = finv(1 - alpha, NumDF, DenDF,0);
    Power = 1 - probf(FCrit, NumDF, DenDF, Noncen);
    call symput('exemplaryPower',power);
  run;
