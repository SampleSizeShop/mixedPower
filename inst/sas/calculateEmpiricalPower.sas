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
%macro longit2Group5Rm(datasetName, covariance);
	proc mixed data=&datasetName;
		model y = trt1_rep1 trt1_rep2 trt1_rep3 trt1_rep4 trt1_rep5
				trt2_rep1 trt2_rep2 trt2_rep3 trt2_rep4 trt2_rep5 / noint solution ddfm=kr;
		repeated / subject=subjectID type=&covariance;
		by setID;
		contrast "time by treatment"
			trt1_rep1 1 trt1_rep2 -1 trt2_rep1 -1 trt2_rep2 1,
			trt1_rep1 1 trt1_rep3 -1 trt2_rep1 -1 trt2_rep3 1,
			trt1_rep1 1 trt1_rep4 -1 trt2_rep1 -1 trt2_rep4 1,
			trt1_rep1 1 trt1_rep5 -1 trt2_rep1 -1 trt2_rep5 1;
	run;
%mend;

* define the mixed model fitting macro for 4 group 5 repeated measures  design; 
%macro longit4Group5Rm(datasetName, covariance);
	proc mixed data=&datasetName;
		model y = trt1_rep1 trt1_rep2 trt1_rep3 trt1_rep4 trt1_rep5
				trt2_rep1 trt2_rep2 trt2_rep3 trt2_rep4 trt2_rep5 
				trt3_rep1 trt3_rep2 trt3_rep3 trt3_rep4 trt3_rep5 
				trt4_rep1 trt4_rep2 trt4_rep3 trt4_rep4 trt4_rep5 
				/ noint solution ddfm=kr;
		repeated / subject=subjectID type=&covariance;
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

* import the parameters defining the designs (generated in R);
proc import datafile="&OUT_DATA_DIR\longitudinalParams.csv"
     out=longitudinalParams
     dbms=csv
     replace;
     getnames=yes;
run;

data longitudinalParams;
	set longitudinalParams(obs=6);
run;

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

	* generate an autoregressive covariance matrix;
	start ar1Matrix(size, rho, matrix);
		dmin=1;
		dmax=size-1;
		ar1 = I(size);
		sizeMinus1 = size-1;
		do r=1 to sizeMinus1;
			rPlus1 = r + 1;
			do c=rPlus1 to size;
				value = rho**(c-r);
		    	ar1[r,c] = value;
		    	ar1[c,r] = value;
			end;
		end;
		matrix = ar1;
	finish;

	* wrapper function to generate appropriate X and Sigma;
	start longitEmpiricalPower(monotone, maxObservations, missingPercent,
								perGroupN, rho, delta, sigmaSq, numGroups, betaScale,
								empiricalPower);

		* Build the sigma matrix for a complete case;
		call ar1Matrix(maxObservations, rho, SigmaComplete);
		SigmaComplete = SigmaComplete # sigmaSq;

		* determine the number of complete and missing cases;
		numComplete = floor(perGroupN * (1 - missingPercent));
		numIncomplete = perGroupN - numComplete;

		* build the design matrix and sigma matrices;
		isu = 1;
		do grp=1 to numgroups;
			designBetween = I(numGroups)[grp,];

			if numIncomplete = 0 then do;
				* add design matrices for each isu;
				X = X // (J(numComplete,1,1) @ (designBetween @ I(maxObservations)));
				* add sigma matrices for each ISU;
				SigmaGroup = I(numComplete) @ SigmaComplete;

				if isu = 1 then SigmaS = SigmaGroup;
				else SigmaS = block(SigmaS, SigmaGroup);

			end;
			else do;
				if monotone = 1 then do;
	      			* monotone dropout pattern - once missing, never come back;
					* 50% complete, 30% missing last observation, 20% missing last 2;
					minObs = maxObservations - 2;
					* add design matrices for each isu;
					X = X // (J(numComplete,1,1) @ (designBetween @ I(maxObservations))) 
						// (J(numIncomplete,1,1) @ (designBetween @ I(maxObservations)[(1:minObs),]));
					* add sigma matrices for each ISU;
					SigmaGroupComplete = I(numComplete) @ SigmaComplete;
					SigmaGroupIncomplete = I(numIncomplete) @ 
						(I(maxObservations)[(1:minObs),] * SigmaComplete * T(I(maxObservations)[(1:minObs),]));
					SigmaGroup = block(SigmaGroupComplete, SigmaGroupIncomplete);

					if isu = 1 then SigmaS = SigmaGroup;
					else SigmaS = block(SigmaS, SigmaGroup);

				end;
				else do;
					* Non-monotone dropout pattern ;
					* - delete 2nd and 4th observations

					* create non-monotone missing patterns;
					missingPattern = {1 3 5}; 
					minObs = 3;
					* add design matrices for each isu;
					X = X // (J(numComplete,1,1) @ (designBetween @ I(maxObservations)))  
						// (J(numIncomplete,1,1) @ (designBetween @ I(maxObservations)[missingPattern,]));
					* add sigma matrices for each ISU;
					SigmaGroupComplete = I(numComplete) @ SigmaComplete;
					SigmaGroupIncomplete = I(numIncomplete) @ (I(maxObservations)[missingPattern,] * SigmaComplete * T(I(maxObservations)[missingPattern,]));
					SigmaGroup = block(SigmaGroupComplete, SigmaGroupIncomplete);
					* concatenate onto full sigma matrix;
					if isu = 1 then SigmaS = SigmaGroup;
					else SigmaS = block(SigmaS, SigmaGroup);

				end;
			end;

			* create subject IDs;
			idList = idList // (T(isu:isu+numComplete-1) @ J(maxObservations,1,1));
			isu = isu + numComplete + 1;
			if numIncomplete > 0 then do;
				idList = idList // (T(isu:isu+numIncomplete-1) @ J(minObs,1,1));
				isu = isu + numIncomplete + 1;
			end;
		end;

		* build column names for X, and build beta;
		if numGroups = 2 then do;
			if maxObservations = 5 then do;
				XFullColNames={"subjectId" "trt1_rep1" "trt1_rep2" "trt1_rep3" "trt1_rep4" "trt1_rep5"
					"trt2_rep1" "trt2_rep2" "trt2_rep3" "trt2_rep4" "trt2_rep5"}; 
				macroName = "longit2Group5Rm";
			end;
			else do;
				XFullColNames={"subjectId" "trt1_rep1" "trt1_rep2" "trt1_rep3"
								"trt2_rep1" "trt2_rep2" "trt2_rep3"}; 
				macroName = "longit2Group3Rm";
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
				XFullColNames={"subjectId" "trt1_rep1" "trt1_rep2" "trt1_rep3" 
				"trt2_rep1" "trt2_rep2" "trt2_rep3" 
				"trt3_rep1" "trt3_rep2" "trt3_rep3" 
				"trt4_rep1" "trt4_rep2" "trt4_rep3" }; 
				macroName = "longit4Group3Rm";
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

		blockSize = 500;
		*if clusterSize > 50 then blockSize = 100;
		call calculateEmpiricalPowerConditional(10000, blockSize,  
		  simlib, simprefix, macroName,
		  X, XFullColNames, XModelColNames, Beta, SigmaS,
		  empiricalPower);


		free X;
		free SigmaS;
	finish;

	* load the parameter data into a matrix;
	paramNames={"monotone" "maxObservations" "missingPercent" "perGroupN" "targetPower" 
				"rho" "delta" "sigmaSq" "numGroups" "betaScale"};
	use longitudinalParams;
  	read all into paramList[colname=paramNames];

	do i=1 to NROW(paramList);
		print ("Case: " + strip(char(i)));

		call longitEmpiricalPower(paramList[i,"monotone"],
							paramList[i,"maxObservations"], 
							paramList[i,"missingPercent"],
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
