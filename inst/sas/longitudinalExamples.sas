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

* define the mixed model fitting macro for 2 group 5 repeated measures design; 
%macro longit2Group5Rm(datasetName);
	proc mixed data=&datasetName;
		model y = trt1_rep1 trt1_rep2 trt1_rep3 trt1_rep4 trt1_rep5
				trt2_rep1 trt2_rep2 trt2_rep3 trt2_rep4 trt2_rep5 / noint solution ddfm=kr;
		repeated / subject=subjectID type=UN;
		by setID;
		contrast "time by treatment"
			trt1_rep1 1 trt1_rep2 -1 trt2_rep1 -1 trt2_rep2 1,
			trt1_rep1 1 trt1_rep3 -1 trt2_rep1 -1 trt2_rep3 1,
			trt1_rep1 1 trt1_rep4 -1 trt2_rep1 -1 trt2_rep4 1,
			trt1_rep1 1 trt1_rep5 -1 trt2_rep1 -1 trt2_rep5 1;
	run;
%mend;

* define the mixed model fitting macro for 2 group 10 repeated measures design; 
%macro longit2Group10Rm(datasetName);
	proc mixed data=&datasetName;
		model y = trt1_rep1 trt1_rep2 trt1_rep3 trt1_rep4 trt1_rep5
				trt1_rep6 trt1_rep7 trt1_rep8 trt1_rep9 trt1_rep10
				trt2_rep1 trt2_rep2 trt2_rep3 trt2_rep4 trt2_rep5
				trt2_rep6 trt2_rep7 trt2_rep8 trt2_rep9 trt2_rep10 / noint solution ddfm=kr;
		repeated / subject=subjectID type=UN;
		by setID;
		contrast "time by treatment"
			trt1_rep1 1 trt1_rep2 -1 trt2_rep1 -1 trt2_rep2 1,
			trt1_rep1 1 trt1_rep3 -1 trt2_rep1 -1 trt2_rep3 1,
			trt1_rep1 1 trt1_rep4 -1 trt2_rep1 -1 trt2_rep4 1,
			trt1_rep1 1 trt1_rep5 -1 trt2_rep1 -1 trt2_rep5 1,
			trt1_rep1 1 trt1_rep6 -1 trt2_rep1 -1 trt2_rep6 1,
			trt1_rep1 1 trt1_rep7 -1 trt2_rep1 -1 trt2_rep7 1,
			trt1_rep1 1 trt1_rep8 -1 trt2_rep1 -1 trt2_rep8 1,
			trt1_rep1 1 trt1_rep9 -1 trt2_rep1 -1 trt2_rep9 1,
			trt1_rep1 1 trt1_rep10 -1 trt2_rep1 -1 trt2_rep10 1;
	run;
%mend;

* define the mixed model fitting macro for 4 group 5 repeated measures  design; 
%macro longit4Group5Rm(datasetName);
	proc mixed data=&datasetName;
		model y = trt1_rep1 trt1_rep2 trt1_rep3 trt1_rep4 trt1_rep5
				trt2_rep1 trt2_rep2 trt2_rep3 trt2_rep4 trt2_rep5 
				trt3_rep1 trt3_rep2 trt3_rep3 trt3_rep4 trt3_rep5 
				trt4_rep1 trt4_rep2 trt4_rep3 trt4_rep4 trt4_rep5 
				/ noint solution ddfm=kr;
		repeated / subject=subjectID type=UN;
		by setID;
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
			trt1_rep1 1 trt1_rep5 -1 trt4_rep1 -1 trt4_rep5 1;
	run;
%mend;

%macro longit4Group10Rm(datasetName);
	proc mixed data=&datasetName;
		model y = trt1_rep1 trt1_rep2 trt1_rep3 trt1_rep4 trt1_rep5
				trt1_rep6 trt1_rep7 trt1_rep8 trt1_rep9 trt1_rep10
				trt2_rep1 trt2_rep2 trt2_rep3 trt2_rep4 trt2_rep5 
				trt2_rep6 trt2_rep7 trt2_rep8 trt2_rep9 trt2_rep10
				trt3_rep1 trt3_rep2 trt3_rep3 trt3_rep4 trt3_rep5 
				trt3_rep6 trt3_rep7 trt3_rep8 trt3_rep9 trt3_rep10
				trt4_rep1 trt4_rep2 trt4_rep3 trt4_rep4 trt4_rep5 
				trt4_rep6 trt4_rep7 trt4_rep8 trt4_rep9 trt4_rep10
				/ noint solution ddfm=kr;
		repeated / subject=subjectID type=UN;
		by setID;
		contrast "time by treatment"
			trt1_rep1 1 trt1_rep2 -1 trt2_rep1 -1 trt2_rep2 1,
			trt1_rep1 1 trt1_rep3 -1 trt2_rep1 -1 trt2_rep3 1,
			trt1_rep1 1 trt1_rep4 -1 trt2_rep1 -1 trt2_rep4 1,
			trt1_rep1 1 trt1_rep5 -1 trt2_rep1 -1 trt2_rep5 1,
			trt1_rep1 1 trt1_rep6 -1 trt2_rep1 -1 trt2_rep6 1,
			trt1_rep1 1 trt1_rep7 -1 trt2_rep1 -1 trt2_rep7 1,
			trt1_rep1 1 trt1_rep8 -1 trt2_rep1 -1 trt2_rep8 1,
			trt1_rep1 1 trt1_rep9 -1 trt2_rep1 -1 trt2_rep9 1,
			trt1_rep1 1 trt1_rep10 -1 trt2_rep1 -1 trt2_rep10 1,

			trt1_rep1 1 trt1_rep2 -1 trt3_rep1 -1 trt3_rep2 1,
			trt1_rep1 1 trt1_rep3 -1 trt3_rep1 -1 trt3_rep3 1,
			trt1_rep1 1 trt1_rep4 -1 trt3_rep1 -1 trt3_rep4 1,
			trt1_rep1 1 trt1_rep5 -1 trt3_rep1 -1 trt3_rep5 1,
			trt1_rep1 1 trt1_rep6 -1 trt3_rep1 -1 trt3_rep6 1,
			trt1_rep1 1 trt1_rep7 -1 trt3_rep1 -1 trt3_rep7 1,
			trt1_rep1 1 trt1_rep8 -1 trt3_rep1 -1 trt3_rep8 1,
			trt1_rep1 1 trt1_rep9 -1 trt3_rep1 -1 trt3_rep9 1,
			trt1_rep1 1 trt1_rep10 -1 trt3_rep1 -1 trt3_rep10 1,

			trt1_rep1 1 trt1_rep2 -1 trt4_rep1 -1 trt4_rep2 1,
			trt1_rep1 1 trt1_rep3 -1 trt4_rep1 -1 trt4_rep3 1,
			trt1_rep1 1 trt1_rep4 -1 trt4_rep1 -1 trt4_rep4 1,
			trt1_rep1 1 trt1_rep5 -1 trt4_rep1 -1 trt4_rep5 1,
			trt1_rep1 1 trt1_rep6 -1 trt4_rep1 -1 trt4_rep6 1,
			trt1_rep1 1 trt1_rep7 -1 trt4_rep1 -1 trt4_rep7 1,
			trt1_rep1 1 trt1_rep8 -1 trt4_rep1 -1 trt4_rep8 1,
			trt1_rep1 1 trt1_rep9 -1 trt4_rep1 -1 trt4_rep9 1,
			trt1_rep1 1 trt1_rep10 -1 trt4_rep1 -1 trt4_rep10 1;
	run;
%mend;


* import the parameters defining the designs (generated in R);
proc import datafile="&OUT_DATA_DIR\longitudinalParams.csv"
     out=longitudinalParams
     dbms=csv
     replace;
     getnames=yes;
run;

data longitudinalParams;
	set longitudinalParams(firstObs=1 obs=1);
run;

proc iml;
maxObservations=5;
min1 = 4;
min2 = 3;
numMin2=1;
numMin1=2;
numComplete=3;
isu = 1;

			* create subject IDs;
			idList = idList // (T(isu:isu+numComplete-1) @ J(maxObservations,1,1));
			isu = isu + numComplete + 1;
			idList = idList // (T(isu:isu+numMin1-1) @ J(min1,1,1));
			isu = isu + numMin1 + 1;
			idList = idList // (T(isu:isu+numMin2-1) @ J(min2,1,1));
			isu = isu + numMin2 + 1;

			print idList;
quit;
/*
* Calculate empirical power for the 4 group, cluster randomized trials
*/
proc iml;
	%INCLUDE "&MODULES_DIR\simulateMixedModel.sxs"/NOSOURCE2;
	%INCLUDE "&MODULES_DIR\calculatePowerKenwardRoger.sxs"/NOSOURCE2;

	* conveience function to generate a lear correlation matrix;
	start learMatrix(size, rho, delta, matrix);
		dmin=1;
		dmax=size-1;
		lear = I(size);
		sizeMinus1 = size-1;
		do r=1 to sizeMinus1;
			rPlus1 = r + 1;
			do c=rPlus1 to size;
				value = rho**(dmin + delta*((c-r-dmin)/(dmax-dmin)));
		    	lear[r,c] = value;
		    	lear[c,r] = value;
			end;
		end;
		matrix = lear;
	finish;

	* wrapper function to generate appropriate X and Sigma;
	start longitEmpiricalPower(monotone, maxObservations, perGroupN, 
							   rho, delta, sigmaSq, numGroups, betaScale,
								empiricalPower);

		* Build the sigma matrix for a complete case;
		call learMatrix(maxObservations, rho, delta, SigmaComplete);
		SigmaComplete = SigmaComplete # sigmaSq;

		* build the design matrix and sigma matrices;
		isu = 1;
		do grp=1 to numgroups;
			designBetween = I(numGroups)[grp,];

			if monotone = 1 then do;
      			* monotone dropout pattern - once missing, never come back;
				* 50% complete, 30% missing last observation, 20% missing last 2;
			    min1 = maxObservations - 1;
				min2 = maxObservations - 2;
				* calculate number of ISUs with each pattern;
				numComplete = floor(0.5 * perGroupN);
				numMin1 = floor(0.3 * perGroupN);
				numMin2 = floor(0.2 * perGroupN);
				* add design matrices for each isu;
				X = X // (J(numComplete,1,1) @ (designBetween @ I(maxObservations))) 
					// (J(numMin1,1,1) @ (designBetween @ I(maxObservations)[(1:min1),])) 
					// (J(numMin2,1,1) @ (designBetween @ I(maxObservations)[(1:min2),]));
				* add sigma matrices for each ISU;
				SigmaGroupComplete = I(numComplete) @ SigmaComplete;
				SigmaGroupMin1 = I(numMin1) @ (I(maxObservations)[(1:min1),] * SigmaComplete * T(I(maxObservations)[(1:min1),]));
				SigmaGroupMin2 = I(numMin2) @ (I(maxObservations)[(1:min2),] * SigmaComplete * T(I(maxObservations)[(1:min2),]));
				SigmaGroup = block(SigmaGroupComplete, SigmaGroupMin1, SigmaGroupMin2);

				if grp = 1 then SigmaS = SigmaGroup;
				else SigmaS = block(SigmaS, SigmaGroup);

			end;
			else do;
				* Non-monotone dropout pattern ;
				* - delete 2nd observations in 30%, delete 3rd in 20%;
				* - if max observations > 3, also delete 4th;

				* calculate number of ISUs with each pattern;
				numComplete = floor(0.5 * perGroupN);
				numMin1 = floor(0.3 * perGroupN);
				numMin2 = floor(0.2 * perGroupN);

				* create non-monotone missing patterns;
				if maxObservations = 3 then do;
					min1Pattern = {1,3};
					min2Pattern = {1,3} || (5:maxObservations);
				end;
				else do;
					min1Pattern = {1,2};
					min2Pattern = {1,2} || (5:maxObservations);
				end;
				
				* add design matrices for each isu;
				X = X // (J(numComplete,1,1) @ (designBetween @ I(maxObservations))) 
					// (J(numMin1,1,1) @ (designBetween @ I(maxObservations)[min1Pattern,])) 
					// (J(numMin2,1,1) @ (designBetween @ I(maxObservations)[min2Pattern,]));
				* add sigma matrices for each ISU;
				SigmaGroupComplete = I(numComplete) @ SigmaComplete;
				SigmaGroupMin1 = I(numMin1) @ (I(maxObservations)[min1Pattern,] * SigmaComplete * T(I(maxObservations)[min1Pattern,]));
				SigmaGroupMin2 = I(numMin2) @ (I(maxObservations)[min2Pattern,] * SigmaComplete * T(I(maxObservations)[min2Pattern,]));
				SigmaGroup = block(SigmaGroupComplete, SigmaGroupMin1, SigmaGroupMin2);
				* concatenate onto full sigma matrix;
				SigmaS = block(SigmaS, SigmaGroup);

			end;

			* create subject IDs;
			idList = idList // (T(isu:isu+numComplete-1) @ J(maxObservations,1,1));
			isu = isu + numComplete + 1;
			idList = idList // (T(isu:isu+numMin1-1) @ J(min1,1,1));
			isu = isu + numMin1 + 1;
			idList = idList // (T(isu:isu+numMin2-1) @ J(min2,1,1));
			isu = isu + numMin2 + 1;
		end;

		* build column names for X, and build beta;
		if numGroups = 2 then do;
			if maxObservations = 5 then do;
				XFullColNames={"subjectId" "trt1_rep1" "trt1_rep2" "trt1_rep3" "trt1_rep4" "trt1_rep5"
					"trt2_rep1" "trt2_rep2" "trt2_rep3" "trt2_rep4" "trt2_rep5"}; 
				macroName = "longit2Group5Rm";
			end;
			else do;
				XFullColNames={"subjectId" "trt1_rep1" "trt1_rep2" "trt1_rep3" "trt1_rep4" "trt1_rep5"
				"trt1_rep6" "trt1_rep7" "trt1_rep8" "trt1_rep9" "trt1_rep10"
				"trt2_rep1" "trt2_rep2" "trt2_rep3" "trt2_rep4" "trt2_rep5"
				"trt2_rep6" "trt2_rep7" "trt2_rep8" "trt2_rep9" "trt2_rep10"}; 
				macroName = "longit2Group10Rm";
			end;
		end;
		else do;
			if maxObservations = 5 then do;
				XFullColNames={"subjectId" "trt1_rep1" "trt1_rep2" "trt1_rep3" "trt1_rep4" "trt1_rep5"
				"trt2_rep1" "trt2_rep2" "trt2_rep3" "trt2_rep4" "trt2_rep5" 
				"trt3_rep1" "trt3_rep2" "trt3_rep3" "trt3_rep4" "trt3_rep5" 
				"trt4_rep1" "trt4_rep2" "trt4_rep3" "trt4_rep4" "trt4_rep5" }; 
				macroName = "longit4Group5Rm";
			end;
			else do;
				XFullColNames={"subjectId" "trt1_rep1" "trt1_rep2" "trt1_rep3" "trt1_rep4" "trt1_rep5"
				"trt1_rep6" "trt1_rep7" "trt1_rep8" "trt1_rep9" "trt1_rep10"
				"trt2_rep1" "trt2_rep2" "trt2_rep3" "trt2_rep4" "trt2_rep5"
				"trt2_rep6" "trt2_rep7" "trt2_rep8" "trt2_rep9" "trt2_rep10"
				"trt3_rep1" "trt3_rep2" "trt3_rep3" "trt3_rep4" "trt3_rep5"
				"trt3_rep6" "trt3_rep7" "trt3_rep8" "trt3_rep9" "trt3_rep10"
				"trt4_rep1" "trt4_rep2" "trt4_rep3" "trt4_rep4" "trt4_rep5" 
				"trt4_rep6" "trt4_rep7" "trt4_rep8" "trt4_rep9" "trt4_rep10"}; 
				macroName = "longit4Group10Rm";
			end;
		end;
		XModelColNames = XFullColNames[1,(2:NCOL(XFullColNames))];
		* define beta;
		Beta = betaScale*({1} // J((NCOL(XModelColNames)-1),1,0));

		print idList;
		X = idList || X ;
		print X;
		print SigmaS;
		
		simlib= "outData";
		simprefix = "longitExamples";

		blockSize = 1000;
		*if clusterSize > 50 then blockSize = 100;
		call calculateEmpiricalPowerConditional(10000, blockSize,  
		  simlib, simprefix, macroName,
		  X, XFullColNames, XModelColNames, Beta, SigmaS,
		  empiricalPower);


		free X;
		free SigmaS;
	finish;

	* load the parameter data into a matrix;
	paramNames={"monotone" "maxObservations" "perGroupN" "targetPower" 
				"rho" "delta" "sigmaSq" "numGroups" "betaScale"};
	use longitudinalParams;
  	read all into paramList[colname=paramNames];

	do i=1 to NROW(paramList);
		print ("Case: " + strip(char(i)));

		call longitEmpiricalPower(paramList[i,"monotone"],
							paramList[i,"maxObservations"], 
							paramList[i,"perGroupN"], 
							paramList[i,"rho"], 
							paramList[i,"delta"], 
							paramList[i,"sigmaSq"],
							paramList[i,"numGroups"],
							paramList[i,"betaScale"],
							empiricalPower);
		empiricalPowerResults = empiricalPowerResults // empiricalPower;
	end;
	
	* create result set;
	resultNames = paramNames || "empiricalPower";
	results = paramList || empiricalPowerResults;
	print results;
	print resultNames;
	* write power results to a data set;
	create longitudinalEmpirical from 
		results[colname=resultNames];
	append from results;
	close longitudinalEmpirical;
quit;

* write the temporary empirical power data set to disk as a csv;
proc export data=longitudinalEmpirical
   outfile="&OUT_DATA_DIR\longitudinalEmpirical.csv"
   dbms=csv
   replace;
run;
