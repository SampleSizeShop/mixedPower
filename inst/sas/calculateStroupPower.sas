/*
*
* Cluster randomized designs with four treatments
*
* We vary the following parameters
*  per group N: 10, 40
*  cluster size: 5, 50, 500
*  missing percent (in 50% of ISUs): 0%, 10%, 20%, 40%
*/

%include "common.sas";


%macro longit2GroupExemplary(betaScale, covariance, missingType, numComplete, numIncomplete);

  * create exemplary data for complete cases;
  data complete;
	do grp = 1 to 2;
	  do isu = 1 to &numComplete; 
	    do rep = 1 to 5;
			y = 3 * (grp = 1 & rep = 1);
			subjectID = ((grp-1) * &numComplete) + isu;
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

  * create exemplary data for incomplete cases;
  %if numIncomplete > 0 %then %do;
	  %if &missingType = "monotone" %then %do;
		  data incomplete;
			do grp = 1 to 2;
			  do isu = (&numComplete * 2) to (&numComplete + &numIncomplete); 
			    do rep = 1 to 3;
					y = 3 * (grp = 1 & rep = 1);
					subjectID = 10 + ((grp-1) * 5) + isu;
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
	  %end;
	  %else %do;
	  	* non-monotone pattern;
	  	  data incomplete;
			do grp = 1 to 2;
			  do isu = 1 to &numComplete; 
			    do rep = 1,3,5;
					y = 3 * (grp = 1 & rep = 1);
					subjectID = 10 + ((grp-1) * 5) + isu;
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
	  %end;

	  * combine the complete and incomplete cases;
	  data exemplaryData;
	    set complete incomplete;
	  run;
  %end;

  * fit the mixed model with the appropriate covariance parameters;
  ods exclude all;
  ods noresults;
  ods output Contrasts=tmpExemplaryContrasts;
  %if &covariance = "CS" %then %do;
    proc glimmix data=exemplaryData;
  		model y = trt1_rep1 trt1_rep2 trt1_rep3 trt1_rep4 trt1_rep5
  				trt2_rep1 trt2_rep2 trt2_rep3 trt2_rep4 trt2_rep5 / noint solution ddfm=kr;
  		random _residual_ / subject=subjectID type=CS;
      parms (1) (0.4) / noiter;
  		contrast "time by treatment"
  			trt1_rep1 1 trt1_rep2 -1 trt2_rep1 -1 trt2_rep2 1,
  			trt1_rep1 1 trt1_rep3 -1 trt2_rep1 -1 trt2_rep3 1,
  			trt1_rep1 1 trt1_rep4 -1 trt2_rep1 -1 trt2_rep4 1,
  			trt1_rep1 1 trt1_rep5 -1 trt2_rep1 -1 trt2_rep5 1;
  	run;
  %end;
  %if &covariance = "CSH" %then %do;
    proc glimmix data=exemplaryData;
    	model y = trt1_rep1 trt1_rep2 trt1_rep3 trt1_rep4 trt1_rep5
  				trt2_rep1 trt2_rep2 trt2_rep3 trt2_rep4 trt2_rep5 / noint solution ddfm=kr;
  		random _residual_ / subject=subjectID type=CSH;
      parms (1) (0.5) (0.3) (0.1) (0.1) (0.4) / noiter;
  		contrast "time by treatment"
  			trt1_rep1 1 trt1_rep2 -1 trt2_rep1 -1 trt2_rep2 1,
  			trt1_rep1 1 trt1_rep3 -1 trt2_rep1 -1 trt2_rep3 1,
  			trt1_rep1 1 trt1_rep4 -1 trt2_rep1 -1 trt2_rep4 1,
  			trt1_rep1 1 trt1_rep5 -1 trt2_rep1 -1 trt2_rep5 1;
  	run;
  %end;
  %if &covariance = "AR(1)" %then %do;
    proc glimmix data=exemplaryData;
    	model y = trt1_rep1 trt1_rep2 trt1_rep3 trt1_rep4 trt1_rep5
  				trt2_rep1 trt2_rep2 trt2_rep3 trt2_rep4 trt2_rep5 / noint solution ddfm=kr;
  		random _residual_ / subject=subjectID type=AR(1);
      parms (1) (0.4) / noiter;
  		contrast "time by treatment"
  			trt1_rep1 1 trt1_rep2 -1 trt2_rep1 -1 trt2_rep2 1,
  			trt1_rep1 1 trt1_rep3 -1 trt2_rep1 -1 trt2_rep3 1,
  			trt1_rep1 1 trt1_rep4 -1 trt2_rep1 -1 trt2_rep4 1,
  			trt1_rep1 1 trt1_rep5 -1 trt2_rep1 -1 trt2_rep5 1;
  	run;
  %end;
  ods results;
  ods exclude none;
  
  * calculate power;
  data power;
    set tmpExemplaryContrasts;
    Noncen = NumDF * FValue;
    Alpha = 0.05;
    FCrit = finv(1 - alpha, NumDF, DenDF,0);
    Power = 1 - probf(FCrit, NumDF, DenDF, Noncen);
    call symput('exemplaryPower',power);
  run;
   
  exemplaryPower;
%mend; /* end longit2GroupExemplary */

%longit2GroupExemplary(3.4, "CS", "monotone", 10, 5);

%macro longit4GroupExemplary(betaScale, covariance, missingType, missingPercent,
                            perGroupN);
  0.15;
%mend; /* end longit4GroupExemplary */

*;
%macro longitExemplary(betaScale, covariance, missingType, missingPercent,
                       perGroupN, numGroups);

	%let numComplete = floor(&perGroupN * (1-&missingPercent));
	%let numIncomplete = &perGroupN - &numComplete;
	%if &numGroups = 2 %then %do; 
		%longit2GroupExemplary(&betaScale, &covariance, &missingType, &numComplete, &numIncomplete);
	%end;
	%if numGroups = 4 %then %do;
		%longit4GroupExemplary(&betaScale, &covariance, &missingType, &numComplete, &numIncomplete);
	%end;
%mend;

data foo;
	exemplaryPower = %longitExemplary(1, "CS", "monotone", 0.2, 50, 2);
run;

/*
* import the parameters defining the designs (generated in R);
proc import datafile="&OUT_DATA_DIR\longitudinalParams.csv"
     out=longitudinalParams
     dbms=csv
     replace;
     getnames=yes;
run;

data longitudinalExemplary;
	set longitudinalParams;
	exemplaryPower = %longitExemplary(betaScale, covariance, missingType, missingPercent,
                       perGroupN, numGroups);
run;


* write the temporary empirical power data set to disk as a csv;
proc export data=longitudinalExemplary
   outfile="&OUT_DATA_DIR\exemplaryPower.csv"
   dbms=csv
   replace;
run;
*/
