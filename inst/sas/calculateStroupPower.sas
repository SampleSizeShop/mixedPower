/*
*
* Calculates power for a set of longitudinal study designs
* using the exemplary data method described by:
* 
* Stroup, W. W. (2013). Generalized linear mixed models: modern concepts, methods and applications.
* (see Chapter 16)
*
* Power results are used in the validation experiment in the manuscript:
*
* Kreidler, S. M., Muller, K. E., & Glueck, D. H. 
* A Power Approximation for Longitudinal Studies Using the 
* Kenward and Roger Wald Test in the Linear Mixed Model, In review.
*/

%include "common.sas";

/*
* Calculate exemplary power for longitudinal designs with 
* 2 treatment groups
*/
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
			  do isu = 1 to &numIncomplete; 
			    do rep = 1 to 3;
					y = &betaScale * (grp = 1 & rep = 1);
					subjectID = (2 * &numComplete) + ((grp-1) * &numIncomplete) + isu;
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
			  do isu = 1 to &numIncomplete; 
			    do rep = 1,3,5;
					y = &betaScale * (grp = 1 & rep = 1);
					subjectID = (2 * &numComplete) + ((grp-1) * &numIncomplete) + isu;
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
	proc glimmix data=exemplaryData noprofile;
		class subjectID;
    	model y = trt1_rep1 trt1_rep2 trt1_rep3 trt1_rep4 trt1_rep5
  				trt2_rep1 trt2_rep2 trt2_rep3 trt2_rep4 trt2_rep5 / noint solution ddfm=bw;
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
    proc glimmix data=exemplaryData noprofile;
		class subjectID;
    	model y = trt1_rep1 trt1_rep2 trt1_rep3 trt1_rep4 trt1_rep5
  				trt2_rep1 trt2_rep2 trt2_rep3 trt2_rep4 trt2_rep5 / noint solution ddfm=bw;
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
    proc glimmix data=exemplaryData noprofile;
		class subjectID;
    	model y = trt1_rep1 trt1_rep2 trt1_rep3 trt1_rep4 trt1_rep5
  				trt2_rep1 trt2_rep2 trt2_rep3 trt2_rep4 trt2_rep5 / noint solution ddfm=bw;
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
    %if %quote(&covariance) = UN %then %do;
    proc glimmix data=exemplaryData noprofile;
		class subjectID;
    	model y = trt1_rep1 trt1_rep2 trt1_rep3 trt1_rep4 trt1_rep5
  				trt2_rep1 trt2_rep2 trt2_rep3 trt2_rep4 trt2_rep5 / noint solution ddfm=bw;
  		random _residual_ / subject=subjectID type=UN;
      parms (1) (0.4) (1) (0.2) (0.7) (1.0) (0.3) (0.6) (0.6) (1.0) (0.2) (0.7) (0.5) (0.8) (1.0)
			/ hold=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15;
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

/*
* Calculate exemplary power for longitudinal designs with 
* 4 treatment groups
*/
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
	  %if %quote(&missingType) = %quote(monotone) %then %do;
		  data incomplete;
			do grp = 1 to 4;
			  do isu = 1 to &numIncomplete; 
			    do rep = 1 to 3;
					y = &betaScale * (grp = 1 & rep = 1);
					subjectID = (4 * &numComplete) + ((grp-1) * &numIncomplete) + isu;
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
			do grp = 1 to 4;
			  do isu = 1 to &numIncomplete; 
			    do rep = 1,3,5;
					y = &betaScale * (grp = 1 & rep = 1);
					subjectID = (4 * &numComplete) + ((grp-1) * &numIncomplete) + isu;
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
	proc glimmix data=exemplaryData noprofile;
		class subjectID;
    	model y = trt1_rep1 trt1_rep2 trt1_rep3 trt1_rep4 trt1_rep5
  				trt2_rep1 trt2_rep2 trt2_rep3 trt2_rep4 trt2_rep5 
				trt3_rep1 trt3_rep2 trt3_rep3 trt3_rep4 trt3_rep5 
				trt4_rep1 trt4_rep2 trt4_rep3 trt4_rep4 trt4_rep5 
		/ noint solution ddfm=bw;
  		random _residual_ / subject=subjectID type=CS;
	    parms (1) (0.4) / hold=1,2;
  		contrast "time by treatment"
  			trt1_rep1 1 trt1_rep2 -1 trt2_rep1 -1 trt2_rep2 1,
  			trt1_rep1 1 trt1_rep3 -1 trt2_rep1 -1 trt2_rep3 1,
  			trt1_rep1 1 trt1_rep4 -1 trt2_rep1 -1 trt2_rep4 1,
  			trt1_rep1 1 trt1_rep5 -1 trt2_rep1 -1 trt2_rep5 1,

			trt1_rep1 1 trt1_rep2 -1 trt3_rep1 -1 trt3_rep2 1,
  			trt1_rep1 1 trt1_rep3 -1 trt3_rep1 -1 trt3_rep3 1,
  			trt1_rep1 1 trt1_rep4 -1 trt3_rep1 -1 trt3_rep4 1,
  			trt1_rep1 1 trt1_rep5 -1 trt3_rep1 -1 trt3_rep5 1,

			trt1_rep1 1 trt1_rep2 -1 trt4_rep1 -1 trt4_rep2 1,
  			trt1_rep1 1 trt1_rep3 -1 trt4_rep1 -1 trt4_rep3 1,
  			trt1_rep1 1 trt1_rep4 -1 trt4_rep1 -1 trt4_rep4 1,
  			trt1_rep1 1 trt1_rep5 -1 trt4_rep1 -1 trt4_rep5 1
		;
	    ods output contrasts=tmpContrasts;
	run;
  %end;
  %if %quote(&covariance) = CSH %then %do;
    proc glimmix data=exemplaryData noprofile;
		class subjectID;
    	model y = trt1_rep1 trt1_rep2 trt1_rep3 trt1_rep4 trt1_rep5
  				trt2_rep1 trt2_rep2 trt2_rep3 trt2_rep4 trt2_rep5 
				trt3_rep1 trt3_rep2 trt3_rep3 trt3_rep4 trt3_rep5 
				trt4_rep1 trt4_rep2 trt4_rep3 trt4_rep4 trt4_rep5 
		/ noint solution ddfm=bw;
  		random _residual_ / subject=subjectID type=CSH;
      parms (1) (0.5) (0.3) (0.1) (0.1) (0.4) / hold=1,2,3,4,5,6;
  		contrast "time by treatment"
  			trt1_rep1 1 trt1_rep2 -1 trt2_rep1 -1 trt2_rep2 1,
  			trt1_rep1 1 trt1_rep3 -1 trt2_rep1 -1 trt2_rep3 1,
  			trt1_rep1 1 trt1_rep4 -1 trt2_rep1 -1 trt2_rep4 1,
  			trt1_rep1 1 trt1_rep5 -1 trt2_rep1 -1 trt2_rep5 1,

			trt1_rep1 1 trt1_rep2 -1 trt3_rep1 -1 trt3_rep2 1,
  			trt1_rep1 1 trt1_rep3 -1 trt3_rep1 -1 trt3_rep3 1,
  			trt1_rep1 1 trt1_rep4 -1 trt3_rep1 -1 trt3_rep4 1,
  			trt1_rep1 1 trt1_rep5 -1 trt3_rep1 -1 trt3_rep5 1,

			trt1_rep1 1 trt1_rep2 -1 trt4_rep1 -1 trt4_rep2 1,
  			trt1_rep1 1 trt1_rep3 -1 trt4_rep1 -1 trt4_rep3 1,
  			trt1_rep1 1 trt1_rep4 -1 trt4_rep1 -1 trt4_rep4 1,
  			trt1_rep1 1 trt1_rep5 -1 trt4_rep1 -1 trt4_rep5 1
		;
		ods output contrasts=tmpContrasts;
  	run;
  %end;
  %if %quote(&covariance) = %quote(AR(1)) %then %do;
    proc glimmix data=exemplaryData noprofile;
		class subjectID;
    	model y = trt1_rep1 trt1_rep2 trt1_rep3 trt1_rep4 trt1_rep5
  				trt2_rep1 trt2_rep2 trt2_rep3 trt2_rep4 trt2_rep5 
				trt3_rep1 trt3_rep2 trt3_rep3 trt3_rep4 trt3_rep5 
				trt4_rep1 trt4_rep2 trt4_rep3 trt4_rep4 trt4_rep5 
		/ noint solution ddfm=bw;
  		random _residual_ / subject=subjectID type=AR(1);
      parms (1) (0.4) / hold=1,2;
  		contrast "time by treatment"
  			trt1_rep1 1 trt1_rep2 -1 trt2_rep1 -1 trt2_rep2 1,
  			trt1_rep1 1 trt1_rep3 -1 trt2_rep1 -1 trt2_rep3 1,
  			trt1_rep1 1 trt1_rep4 -1 trt2_rep1 -1 trt2_rep4 1,
  			trt1_rep1 1 trt1_rep5 -1 trt2_rep1 -1 trt2_rep5 1,

			trt1_rep1 1 trt1_rep2 -1 trt3_rep1 -1 trt3_rep2 1,
  			trt1_rep1 1 trt1_rep3 -1 trt3_rep1 -1 trt3_rep3 1,
  			trt1_rep1 1 trt1_rep4 -1 trt3_rep1 -1 trt3_rep4 1,
  			trt1_rep1 1 trt1_rep5 -1 trt3_rep1 -1 trt3_rep5 1,

			trt1_rep1 1 trt1_rep2 -1 trt4_rep1 -1 trt4_rep2 1,
  			trt1_rep1 1 trt1_rep3 -1 trt4_rep1 -1 trt4_rep3 1,
  			trt1_rep1 1 trt1_rep4 -1 trt4_rep1 -1 trt4_rep4 1,
  			trt1_rep1 1 trt1_rep5 -1 trt4_rep1 -1 trt4_rep5 1
		;
		ods output contrasts=tmpContrasts;
  	run;
  %end;
    %if %quote(&covariance) = UN %then %do;
    proc glimmix data=exemplaryData noprofile;
		class subjectID;
    	model y = trt1_rep1 trt1_rep2 trt1_rep3 trt1_rep4 trt1_rep5
  				trt2_rep1 trt2_rep2 trt2_rep3 trt2_rep4 trt2_rep5 
				trt3_rep1 trt3_rep2 trt3_rep3 trt3_rep4 trt3_rep5 
				trt4_rep1 trt4_rep2 trt4_rep3 trt4_rep4 trt4_rep5 
		/ noint solution ddfm=bw;
  		random _residual_ / subject=subjectID type=UN;
      parms (1) (0.4) (1) (0.2) (0.7) (1.0) (0.3) (0.6) (0.6) (1.0) (0.2) (0.7) (0.5) (0.8) (1.0)
			/ hold=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15;
  		contrast "time by treatment"
  			trt1_rep1 1 trt1_rep2 -1 trt2_rep1 -1 trt2_rep2 1,
  			trt1_rep1 1 trt1_rep3 -1 trt2_rep1 -1 trt2_rep3 1,
  			trt1_rep1 1 trt1_rep4 -1 trt2_rep1 -1 trt2_rep4 1,
  			trt1_rep1 1 trt1_rep5 -1 trt2_rep1 -1 trt2_rep5 1,

			trt1_rep1 1 trt1_rep2 -1 trt3_rep1 -1 trt3_rep2 1,
  			trt1_rep1 1 trt1_rep3 -1 trt3_rep1 -1 trt3_rep3 1,
  			trt1_rep1 1 trt1_rep4 -1 trt3_rep1 -1 trt3_rep4 1,
  			trt1_rep1 1 trt1_rep5 -1 trt3_rep1 -1 trt3_rep5 1,

			trt1_rep1 1 trt1_rep2 -1 trt4_rep1 -1 trt4_rep2 1,
  			trt1_rep1 1 trt1_rep3 -1 trt4_rep1 -1 trt4_rep3 1,
  			trt1_rep1 1 trt1_rep4 -1 trt4_rep1 -1 trt4_rep4 1,
  			trt1_rep1 1 trt1_rep5 -1 trt4_rep1 -1 trt4_rep5 1
		;
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
	*delete complete incomplete;
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
	set longitudinalParams;*(firstobs=136 obs=136);
	call execute('%longitExemplary(' || _N_ || ',' || betaScale || ',' || covariance || ',' || missingType || 
					',' || missingPercent || ',' || perGroupN || ',' || numGroups || ', exemplaryPower)');
run;

* merge the power results with the input parameters;
data longitudinalExemplaryPower;
	set longitudinalParams;
	set exemplaryPower(keep=exemplaryPower);	

run;
/*
proc freq data=longitudinalParams;
	table covariance;
run;
*/

* write the temporary empirical power data set to disk as a csv;
proc export data=longitudinalExemplaryPower
   outfile="&OUT_DATA_DIR\longitudinalExemplaryPower.csv"
   dbms=csv
   replace;
run;

