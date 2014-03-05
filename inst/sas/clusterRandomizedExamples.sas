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

* define the mixed model fitting macro for 4 group design; 
%macro clustered4Group(datasetName);
	proc mixed data=&datasetName;
		model y = trt1 trt2 trt3 trt4 / noint solution;
		random int / subject=clusterID;
		by setID;
		contrast "4 group main effect" 
			trt1 1 trt2 -1,
			trt1 1 trt3 -1,
			trt1 1 trt4 -1;
	run;
%mend;

* define the mixed model fitting macro for 2 group design; 
%macro clustered2Group(datasetName);
	proc mixed data=&datasetName;
		model y = trt1 trt2 / noint solution;
		random int / subject=clusterID;
		by setID;
		contrast "2 group main effect" trt1 1 trt2 -1;
	run;
%mend;

* import the parameters defining the designs (generated in R);
proc import datafile="&OUT_DATA_DIR\clusterRandomizedParams.csv"
     out=clusterRandomizedParams
     dbms=csv
     replace;
     getnames=yes;
run;


data clusterRandomizedParams;
	set clusterRandomizedParams(obs=3);
run;

/*
* Calculate empirical power for the 4 group, cluster randomized trials
*/
proc iml;
	%INCLUDE "&MODULES_DIR\simulateMixedModel.sxs"/NOSOURCE2;
	%INCLUDE "&MODULES_DIR\calculatePowerKenwardRoger.sxs"/NOSOURCE2;

	* wrapper function to generate appropriate X and Sigma;
	start clusterEmpiricalPower(numGroups, perGroupN, clusterSize, 
								deletePercent, sigmaSq, icc,
								empiricalPower);
		* cell means essence matrix;
		Xessence = I(numGroups);

		if numGroups = 2 then do;
			XFullColNames={"clusterId" "trt1" "trt2"}; 
			XModelColNames = {"trt1" "trt2"};
			* define beta;
			Beta = {1,0};
			macroName = "clustered2Group";
		end;
		else do;
			XFullColNames={"clusterId" "trt1" "trt2" "trt3" "trt4"}; 
			XModelColNames = {"trt1" "trt2" "trt3" "trt4"};
			* define beta;
			Beta = {1,0,0,0};
			macroName = "clustered4Group";
		end;

		* calculate the number of complete cases;
		numComplete = perGroupN / 2;
		numIncomplete = perGroupN - numComplete;
		* calculate the size of complete and incomplete cases;
		incompleteSize = floor(clusterSize*(1-deletePercent));
		* build the X matrix;
		totalObsPerGroup = numComplete*clusterSize + numIncomplete*incompleteSize;
		X = Xessence@J(totalObsPerGroup,1,1);
		* build the cluster id's;
		do grp=1 to numGroups;
			do isu=1 to perGroupN by 2;
			    id = (grp-1)*perGroupN + isu;
				idList = idList // J(clusterSize,1,id) //
					J(incompleteSize,1,id + 1);
			end;
		end;
		X = idList || X ;
		*print X;

		SigmaComplete = sigmaSq#(J(clusterSize,clusterSize,1)*icc+I(clusterSize)*(1-icc));
		SigmaIncomplete = sigmaSq#(J(incompleteSize,incompleteSize,1)*icc+I(incompleteSize)*(1-icc));
		SigmaS = I((perGroupN/2)*numGroups) @ block(SigmaComplete, SigmaIncomplete);
		*print SigmaS;
		
		simlib= "outData";
		simprefix = "clusterExamples";

		call calculateEmpiricalPowerConditional(10000, 1000,  
		  simlib, simprefix, macroName,
		  X, XFullColNames, XModelColNames, Beta, SigmaS,
		  empiricalPower);


		free X;
		free SigmaS;
	finish;

	* load the parameter data into a matrix;
	paramNames={"missingPercent","clusterSize","perGroupN","targetPower",
				"sigmaSq","icc","numGroups","betaScale"};
	use clusterRandomizedParams;
  	read all into paramList[colname=paramNames];
	
	do i=1 to NROW(paramList);
		print ("Case: " + strip(char(i)));
		call clusterEmpiricalPower(paramList[i,"numGroups"],
							paramList[i,"perGroupN"], 
							paramList[i,"clusterSize"], 
							paramList[i,"missingPercent"], 
							paramList[i,"sigmaSq"], 
							paramList[i,"icc"],
							empiricalPower);
		empiricalPowerResults = empiricalPowerResults // empiricalPower;
	end;
	
	* create result set;
	resultNames = paramNames || "empiricalPower";
	results = paramList || empiricalPowerResults;
	print results;
	print resultNames;
	* write power results to a data set;
	create clusterRandomizedEmpirical from 
		results[colname=resultNames];
	append from results;
	close clusterRandomizedEmpirical;
quit;

* write the temporary empirical power data set to disk as a csv;
proc export data=clusterFourGroupEmpirical
   outfile='&OUTDATA\clusterFourGroupEmpirical.csv'
   dbms=csv
   replace;
run;
