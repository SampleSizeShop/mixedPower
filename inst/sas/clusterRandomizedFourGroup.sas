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

* define the mixed model fitting macro; 
* this must contain a 'by setID' statement, but can otherwise;
* be defined as needed by the model;
* define the mixed model fitting macro; 
* this must contain a 'by setID' statement, but can otherwise;
* be defined as needed by the model;
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

/*
* Calculate empirical power for the 4 group, cluster randomized trials
*/
proc iml;
	%INCLUDE "&MODULES_DIR\simulateMixedModel.sxs"/NOSOURCE2;
	%INCLUDE "&MODULES_DIR\calculatePowerKenwardRoger.sxs"/NOSOURCE2;

	* wrapper function to generate appropriate X and Sigma;
	start clusterEmpiricalPower(perGroupN, clusterSize, deletePercent, empiricalPower);
		* cell means essence matrix;
		Xessence = I(4);
		XFullColNames={"clusterId" "trt1" "trt2" "trt3" "trt4"}; 
		XModelColNames = {"trt1" "trt2" "trt3" "trt4"};
		* calculate the number of complete cases;
		numComplete = perGroupN / 2;
		numIncomplete = perGroupN - numComplete;
		* calculate the size of complete and incomplete cases;
		incompleteSize = floor(clusterSize*(1-deletePercent));
		* build the X matrix;
		totalObsPerGroup = numComplete*clusterSize + numIncomplete*incompleteSize;
		X = Xessence@J(totalObsPerGroup,1,1);
		* build the cluster id's;
		do grp=1 to 4;
			do isu=1 to perGroupN by 2;
			    id = (grp-1)*perGroupN + isu;
				idList = idList // J(clusterSize,1,id) //
					J(incompleteSize,1,id + 1);
			end;
		end;
		X = idList || X ;
		*print X;

		* define beta;
		Beta = {1,0,0,0};
		SigmaComplete = 2#(J(clusterSize,clusterSize,1)*0.04+I(clusterSize)*0.96);
		SigmaIncomplete = 2#(J(incompleteSize,incompleteSize,1)*0.04+I(incompleteSize)*0.96);
		SigmaS = I((perGroupN/2)*4) @ block(SigmaComplete, SigmaIncomplete);
		*print SigmaS;
		
		simlib= "outData";
		simprefix = "crT4N" + strip(char(perGroupN)) + 
				"C" + strip(char(clusterSize)) + "D" + strip(char(deletePercent));

		call calculateEmpiricalPowerConditional(10000, 1000,  
		  simlib, simprefix, "clustered4Group",
		  X, XFullColNames, XModelColNames, Beta, SigmaS,
		  empiricalPower);


		free X;
		free SigmaS;
	finish;

	* parameters to vary;
	* number of clusters in each treatment group;
	perGroupNList = {10};* 40};
	* cluster size for complete case;
	clusterSizeList = {5};* 50 100 500};
	* delete percent for incomplete cases;
	deletePercentList = {0};* 0.1 0.2 0.4};

	* run the cases ;
	caseCounter = 1;
	do perGroupNIdx=1 to NCOL(perGroupNList);
		do clusterSizeIdx=1 to NCOL(clusterSizeList);
			do deletePercentIdx=1 to NCOL(deletePercentList);
				print ("Case: " + strip(char(caseCounter)));
				call clusterEmpiricalPower(perGroupNList[perGroupNIdx], 
									clusterSizeList[clusterSizeIdx], 
									deletePercentList[deletePercentIdx],
									empiricalPower);
									print(empiricalPower);
				empiricalPowerResults = empiricalPowerResults // 
					(perGroupN || clusterSize || deletePercent || empiricalPower);
				caseCounter = caseCounter + 1;
			end;
		end;
	end;
	
	* write power results to a data set;
	create clusterFourGroupEmpirical from 
		empiricalPowerResults[colname={"perGroupN" "clusterSize" "deletePercent" "empiricalPower"}];
	append from empiricalPowerResults;
	close clusterFourGroupEmpirical;
quit;

* write the temporary empirical power data set to disk as a csv;
proc export data=clusterFourGroupEmpirical
   outfile='&OUTDATA\clusterFourGroupEmpirical.csv'
   dbms=csv
   replace;
run;
