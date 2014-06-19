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


%macro longit2GroupExemplary(case, betaScale, covariance, missingType, numComplete, numIncomplete,
	powerDataSet);

  * create exemplary data for complete cases;
  data complete;
	do grp = 1 to 2;
	  do isu = 1 to &numComplete; 
	    do rep = 1 to 5;
			y = &betaScale * (grp = 1 & rep = 1);
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

  data exemplaryData;
    set complete;
  run;

  * create exemplary data for incomplete cases;
  %if &numIncomplete > 0 %then %do;
	  %if %quote(&missingType) = %quote(monotone) %then %do;
		  data incomplete;
			do grp = 1 to 2;
			  do isu = (&numComplete * 2) to ((&numComplete*2) + &numIncomplete); 
			    do rep = 1 to 3;
					y = &betaScale * (grp = 1 & rep = 1);
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
			  do isu = (&numComplete * 2) to ((&numComplete*2) + &numIncomplete); 
			    do rep = 1,3,5;
					y = &betaScale * (grp = 1 & rep = 1);
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
  %if %quote(&covariance) = CS %then %do;
	proc glimmix data=exemplaryData;
		class subjectID;
    	model y = trt1_rep1 trt1_rep2 trt1_rep3 trt1_rep4 trt1_rep5
  				trt2_rep1 trt2_rep2 trt2_rep3 trt2_rep4 trt2_rep5 / noint solution ddfm=residual;
  		random _residual_ / subject=subjectID type=CS;
	    parms (1) (0.4) / hold=1,2;
  		contrast "time by treatment"
  			trt1_rep1 1 trt1_rep2 -1 trt2_rep1 -1 trt2_rep2 1,
  			trt1_rep1 1 trt1_rep3 -1 trt2_rep1 -1 trt2_rep3 1,
  			trt1_rep1 1 trt1_rep4 -1 trt2_rep1 -1 trt2_rep4 1,
  			trt1_rep1 1 trt1_rep5 -1 trt2_rep1 -1 trt2_rep5 1;
	    ods output contrasts=tmpContrasts;
	run;
  %end;
  %if %quote(&covariance) = CSH %then %do;
    proc glimmix data=exemplaryData;
		class subjectID;
    	model y = trt1_rep1 trt1_rep2 trt1_rep3 trt1_rep4 trt1_rep5
  				trt2_rep1 trt2_rep2 trt2_rep3 trt2_rep4 trt2_rep5 / noint solution ddfm=residual;
  		random _residual_ / subject=subjectID type=CSH;
      parms (1) (0.5) (0.3) (0.1) (0.1) (0.4) / hold=1,2,3,4,5,6;
  		contrast "time by treatment"
  			trt1_rep1 1 trt1_rep2 -1 trt2_rep1 -1 trt2_rep2 1,
  			trt1_rep1 1 trt1_rep3 -1 trt2_rep1 -1 trt2_rep3 1,
  			trt1_rep1 1 trt1_rep4 -1 trt2_rep1 -1 trt2_rep4 1,
  			trt1_rep1 1 trt1_rep5 -1 trt2_rep1 -1 trt2_rep5 1;
		ods output contrasts=tmpContrasts;
  	run;
  %end;
  %if %quote(&covariance) = %quote(AR(1)) %then %do;
    proc glimmix data=exemplaryData;
		class subjectID;
    	model y = trt1_rep1 trt1_rep2 trt1_rep3 trt1_rep4 trt1_rep5
  				trt2_rep1 trt2_rep2 trt2_rep3 trt2_rep4 trt2_rep5 / noint solution ddfm=residual;
  		random _residual_ / subject=subjectID type=AR(1);
      parms (1) (0.4) / hold=1,2;
  		contrast "time by treatment"
  			trt1_rep1 1 trt1_rep2 -1 trt2_rep1 -1 trt2_rep2 1,
  			trt1_rep1 1 trt1_rep3 -1 trt2_rep1 -1 trt2_rep3 1,
  			trt1_rep1 1 trt1_rep4 -1 trt2_rep1 -1 trt2_rep4 1,
  			trt1_rep1 1 trt1_rep5 -1 trt2_rep1 -1 trt2_rep5 1;
		ods output contrasts=tmpContrasts;
  	run;
  %end;
  ods results;
  ods exclude none;
  
  * calculate power;
  data power;
    set tmpContrasts;
	case = &case;
    Noncen = NumDF * FValue;
    Alpha = 0.05;
    FCrit = finv(1 - alpha, NumDF, DenDF,0);
    exemplaryPower = 1 - probf(FCrit, NumDF, DenDF, Noncen);
  run;

  PROC datasets; 
  	append base=&powerDataSet data=power force; 
	delete complete incomplete;
  quit;
   
%mend; /* end longit2GroupExemplary */

%macro longit4GroupExemplary(case, betaScale, covariance, missingType, missingPercent,
                            perGroupN, powerDataSet);
	%put complete: &numComplete, incomplete: &numIncomplete;

  proc datasets;
	delete complete incomplete;
  quit;
  * create exemplary data for complete cases;
  data complete;
	do grp = 1 to 4;
	  do isu = 1 to &numComplete; 
	    do rep = 1 to 5;
			y = &betaScale * (grp = 1 & rep = 1);
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
			trt3_rep1 = (grp = 3 & rep = 1); 
			trt3_rep2 = (grp = 3 & rep = 2);
			trt3_rep3 = (grp = 3 & rep = 3);
			trt3_rep4 = (grp = 3 & rep = 4);
			trt3_rep5 = (grp = 3 & rep = 5);
			trt4_rep1 = (grp = 4 & rep = 1); 
			trt4_rep2 = (grp = 4 & rep = 2);
			trt4_rep3 = (grp = 4 & rep = 3);
			trt4_rep4 = (grp = 4 & rep = 4);
			trt4_rep5 = (grp = 4 & rep = 5);
			output;
		end;
	  end;
    end;
	drop isu grp rep;
  run;

  data exemplaryData;
    set complete;
  run;

  * create exemplary data for incomplete cases;
  %if &numIncomplete > 0 %then %do;
	  %if &missingType = "monotone" %then %do;
		  data incomplete;
			do grp = 1 to 2;
			  do isu = (&numComplete * 4) to ((&numComplete*4) + &numIncomplete); 
			    do rep = 1 to 3;
					y = &betaScale * (grp = 1 & rep = 1);
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
					trt3_rep1 = (grp = 3 & rep = 1); 
					trt3_rep2 = (grp = 3 & rep = 2);
					trt3_rep3 = (grp = 3 & rep = 3);
					trt3_rep4 = (grp = 3 & rep = 4);
					trt3_rep5 = (grp = 3 & rep = 5);
					trt4_rep1 = (grp = 4 & rep = 1); 
					trt4_rep2 = (grp = 4 & rep = 2);
					trt4_rep3 = (grp = 4 & rep = 3);
					trt4_rep4 = (grp = 4 & rep = 4);
					trt4_rep5 = (grp = 4 & rep = 5);
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
			  do isu = (&numComplete * 4) to ((&numComplete*4) + &numIncomplete); 
			    do rep = 1,3,5;
					y = &betaScale * (grp = 1 & rep = 1);
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
					trt3_rep1 = (grp = 3 & rep = 1); 
					trt3_rep2 = (grp = 3 & rep = 2);
					trt3_rep3 = (grp = 3 & rep = 3);
					trt3_rep4 = (grp = 3 & rep = 4);
					trt3_rep5 = (grp = 3 & rep = 5);
					trt4_rep1 = (grp = 4 & rep = 1); 
					trt4_rep2 = (grp = 4 & rep = 2);
					trt4_rep3 = (grp = 4 & rep = 3);
					trt4_rep4 = (grp = 4 & rep = 4);
					trt4_rep5 = (grp = 4 & rep = 5);
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
  %if %quote(&covariance) = CS %then %do;
	proc glimmix data=exemplaryData;
		class subjectID;
    	model y = trt1_rep1 trt1_rep2 trt1_rep3 trt1_rep4 trt1_rep5
  				trt2_rep1 trt2_rep2 trt2_rep3 trt2_rep4 trt2_rep5 / noint solution ddfm=residual;
  		random _residual_ / subject=subjectID type=CS;
	    parms (1) (0.4) / hold=1,2;
  		contrast "time by treatment"
  			trt1_rep1 1 trt1_rep2 -1 trt2_rep1 -1 trt2_rep2 1,
  			trt1_rep1 1 trt1_rep3 -1 trt2_rep1 -1 trt2_rep3 1,
  			trt1_rep1 1 trt1_rep4 -1 trt2_rep1 -1 trt2_rep4 1,
  			trt1_rep1 1 trt1_rep5 -1 trt2_rep1 -1 trt2_rep5 1;
	    ods output contrasts=tmpContrasts;
	run;
  %end;
  %if %quote(&covariance) = CSH %then %do;
    proc glimmix data=exemplaryData;
		class subjectID;
    	model y = trt1_rep1 trt1_rep2 trt1_rep3 trt1_rep4 trt1_rep5
  				trt2_rep1 trt2_rep2 trt2_rep3 trt2_rep4 trt2_rep5 / noint solution ddfm=residual;
  		random _residual_ / subject=subjectID type=CSH;
      parms (1) (0.5) (0.3) (0.1) (0.1) (0.4) / hold=1,2,3,4,5,6;
  		contrast "time by treatment"
  			trt1_rep1 1 trt1_rep2 -1 trt2_rep1 -1 trt2_rep2 1,
  			trt1_rep1 1 trt1_rep3 -1 trt2_rep1 -1 trt2_rep3 1,
  			trt1_rep1 1 trt1_rep4 -1 trt2_rep1 -1 trt2_rep4 1,
  			trt1_rep1 1 trt1_rep5 -1 trt2_rep1 -1 trt2_rep5 1;
		ods output contrasts=tmpContrasts;
  	run;
  %end;
  %if %quote(&covariance) = %quote(AR(1)) %then %do;
    proc glimmix data=exemplaryData;
		class subjectID;
    	model y = trt1_rep1 trt1_rep2 trt1_rep3 trt1_rep4 trt1_rep5
  				trt2_rep1 trt2_rep2 trt2_rep3 trt2_rep4 trt2_rep5 / noint solution ddfm=residual;
  		random _residual_ / subject=subjectID type=AR(1);
      parms (1) (0.4) / hold=1,2;
  		contrast "time by treatment"
  			trt1_rep1 1 trt1_rep2 -1 trt2_rep1 -1 trt2_rep2 1,
  			trt1_rep1 1 trt1_rep3 -1 trt2_rep1 -1 trt2_rep3 1,
  			trt1_rep1 1 trt1_rep4 -1 trt2_rep1 -1 trt2_rep4 1,
  			trt1_rep1 1 trt1_rep5 -1 trt2_rep1 -1 trt2_rep5 1;
		ods output contrasts=tmpContrasts;
  	run;
  %end;
  ods results;
  ods exclude none;
  
  * calculate power;
  data power;
    set tmpContrasts;
	case = &case;
    Noncen = NumDF * FValue;
    Alpha = 0.05;
    FCrit = finv(1 - alpha, NumDF, DenDF,0);
    exemplaryPower = 1 - probf(FCrit, NumDF, DenDF, Noncen);
  run;

  PROC datasets; 
  	append base=&powerDataSet data=power force; 
	delete complete incomplete;
  quit;

%mend; /* end longit4GroupExemplary */

*;
%macro longitExemplary(case, betaScale, covariance, missingType, missingPercent,
                       perGroupN, numGroups, powerDataSet);

	%let numComplete = %sysevalf(&perGroupN * (1-&missingPercent), floor);
	%let numIncomplete = %sysevalf(&perGroupN - &numComplete,integer);

%put Complete: &numComplete &numIncomplete;
		
	%if &numGroups = 2 %then %do; 
		%longit2GroupExemplary(&case, &betaScale, &covariance, &missingType, &numComplete, &numIncomplete,
			&powerDataSet);
	%end;
	%if &numGroups = 4 %then %do;
		%longit4GroupExemplary(&case, &betaScale, &covariance, &missingType, &numComplete, &numIncomplete,
			&powerDataSet);
	%end;
%mend;

* import the parameters defining the designs (generated in R);
proc import datafile="&OUT_DATA_DIR\longitudinalParams.csv"
     out=longitudinalParams
     dbms=csv
     replace;
     getnames=yes;
run;

proc datasets;
	delete exemplaryPower;
quit;
data _null_;
	set longitudinalParams;*(firstobs=109 obs=109);
	call execute('%longitExemplary(' || _N_ || ',' || betaScale || ',' || covariance || ',' || missingType || 
					',' || missingPercent || ',' || perGroupN || ',' || numGroups || ', exemplaryPower)');
run;

/*
proc freq data=longitudinalParams;
	table covariance;
run;
*/
/*
* write the temporary empirical power data set to disk as a csv;
proc export data=longitudinalExemplary
   outfile="&OUT_DATA_DIR\exemplaryPower.csv"
   dbms=csv
   replace;
run;
*/
