/*
*
* Calculate empirical power values for the validation study in
* the manuscript in README
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
	set longitudinalParams;
	* make indicators for character data;
	monotone = (missingType = "monotone");
	covarCS = (covariance = "CS");
	covarCSH = (covariance = "CSH");
	covarAR1 = (covariance = "AR(1)");
	covarUN = (covariance = "UN");
run;

/*
* Calculate empirical power for the 4 group, cluster randomized trials
*/
proc iml;
	* conveience function to generate a lear correlation matrix;
	start learMatrix(size, rho, delta);
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
		return(lear);
	finish;

	* generate an autoregressive covariance matrix;
	start ar1Matrix(size, rho);
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
		return(ar1);
	finish;

	/* 
	* generate a heterogeneous compound symmetric covariance
	* sigmaList - an Nx1 vector of sigma values
	* rho - the correlation between observations
	*/
	start heterogeneousCSMatrix(sigmaList,rho);
		size = nrow(sigmaList);
		matrix = J(size,size,0);
		do r = 1 to size;
			do c = 1 to size;
		      if r = c then
		        matrix[r,c] = sigmaList[r]*sigmaList[c];
		      else
		        matrix[r,c] = sigmaList[r]*sigmaList[c] * rho;
			end;
		end;
		return(matrix);
	finish;


	/*
	* Generate data for a longitudinal design with the specified parameters
	*/
	start generateData(replicates, setSize, simLib, simPrefix,
			monotone, covarCS, covarCSH, covarAR1,
			maxObservations, missingPercent,
			perGroupN, numGroups, betaScale, dsList);

		* determine the number of complete and incomplete sampling units;
		numComplete = floor(perGroupN * (1 - missingPercent));
		numIncomplete = perGroupN - numComplete;

		* build the within design matrix for complete cases and incomplete cases;
		designWithinComplete = I(maxObservations);
		X = J(numComplete,1,1) @ I(numGroups) @ designWithinComplete; 
		if numIncomplete > 0 then do;
			if monotone = 1 then
			    designWithinIncomplete = I(maxObservations)[1:(maxObservations-2),];
			else
				designWithinIncomplete = I(maxObservations)[{1,3,5},];
			X = X // J(numIncomplete,1,1) @ I(numGroups) @ designWithinIncomplete; 
		end;

		* build the beta matrix;
		Beta = {1} // J(((numGroups * maxObservations) - 1), 1, 0);
		Beta = Beta * betaScale;

		* create the Sigma matrix for the complete cases and incomplete cases; 
		rho = 0.04;
		if covarCS then
			SigmaComplete = heterogeneousCSMatrix(J(maxObservations,1,1), rho);
		else do;
			if covarCSH then
				SigmaComplete = heterogeneousCSMatrix({1.0,0.5,0.3,0.1,0.1}, rho);
			else do;
				if covarAR1 then
					SigmaComplete = ar1Matrix(maxObservations, rho);
				else
					SigmaComplete = {
						1 0.4 0.2 0.3 0.2,
                        0.4 1 0.7 0.6 0.7,
                        0.2 0.7 1 0.6 0.5,
                        0.3 0.6 0.6 1 0.8,
                        0.2 0.7 0.5 0.8 1
						};
			end;
		end;
		if numIncomplete > 0 then
			SigmaIncomplete = designWithinIncomplete * SigmaComplete * designWithinIncomplete`;

			*print SigmaIncomplete;
		* build the data set column names;
		if numGroups = 2 then
			XFullColNames={"subjectId" "trt1_rep1" "trt1_rep2" "trt1_rep3" "trt1_rep4" "trt1_rep5"
					"trt2_rep1" "trt2_rep2" "trt2_rep3" "trt2_rep4" "trt2_rep5"}; 
		else
			XFullColNames={"subjectId" "trt1_rep1" "trt1_rep2" "trt1_rep3" "trt1_rep4" "trt1_rep5"
			"trt2_rep1" "trt2_rep2" "trt2_rep3" "trt2_rep4" "trt2_rep5" 
			"trt3_rep1" "trt3_rep2" "trt3_rep3" "trt3_rep4" "trt3_rep5" 
			"trt4_rep1" "trt4_rep2" "trt4_rep3" "trt4_rep4" "trt4_rep5" }; 

		* determine the number of sets;
		numSets = replicates / setSize;

		* create a vector of means (all 0's) for the random errors;
		mu = J(maxObservations,1,0);
		if numIncomplete > 0 then do;
			if monotone = 1 then
			    muIncomplete = J(maxObservations-2,1,0);
			else
				muIncomplete = J(3,1,0);
		end;

		/*
		*
		* We create a total of #replicates data sets in the
		* long format appropriate for PROC MIXED.
		*
		* The replicates are split across multiple files of size 'setSize'
		* with the start and end replicate number in the filename.
		*
		* Each data contains the following fields
		*  setID - the data set identifier
		*  Y - the simulated Y values
		*  E - the simulated errors
		*  X - all columns from the user specified X matrix with
		*    column names as specified in XFullColNames
		*/
		setId = 1;		

		do setNum = 1 to numSets;
		  do repNum = 1 to setSize;
		    isu = 1;
		    * generate complete cases;
		    E = E // colvec(randnormal(numComplete * numGroups, mu, SigmaComplete));
			idList = idList // (T(isu:isu+(numComplete*numGroups)-1) @ J(maxObservations,1,1));
			isu = isu + (numComplete*numGroups);

			if numIncomplete > 0 then do;
			    E = E // colvec(randnormal(numIncomplete * numGroups, muIncomplete, SigmaIncomplete));
				idList = idList // (T(isu:isu+(numIncomplete*numGroups)-1) @ J(nrow(SigmaIncomplete),1,1));
			end;
			* form the responses;
			Y = X * Beta + E;

			* build the data set;
			block = J(nrow(Y),1,setId) || Y || E || idList || X;
			setId = setId + 1;

			* append the block to the data set;
			output = output // block;

			free E;
			free Y;
			free idList;
		  end;
		  * write to disk;
		  startIterNum = (setNum-1)*setSize+1;
		  startIter = char(startIterNum);
		  endIterNum = startIterNum + setSize - 1;
		  endIter = char(endIterNum); 
		  dataSetName = simPrefix + "Iter" + strip(startIter) + 
				"to" + strip(endIter);

		  dsList = dsList // dataSetname;

		  /*
		  * Write a temporary data set since I
		  * can't convince IML to use a dynamic name here
		  */
		  names = ({"setID" "Y" "E"}) || XFullColNames;
		  create temp from output[colname=names];
		  append from output;
		  close temp;
		  free output;

		  /*
		  * Change the name of the data set using SAS
		  */
		  submit dataSetName;
			data &dataSetName;
				set temp;
			run;
		  endsubmit;
		  
		  * remove the temp data set;
		  call delete("work", "temp"); 
		end;

	finish;


	* wrapper function to generate appropriate X and Sigma;
	start longitEmpiricalPower(replicates, setSize, simlib, simprefix,
			monotone, covarCS, covarCSH, covarAR1, 
			maxObservations, missingPercent,
			perGroupN, numGroups, betaScale, empiricalPower);

		* generate data for the design;
		run generateData(replicates, setSize, simlib, simprefix,
			monotone, covarCS, covarCSH, covarAR1,
			maxObservations, missingPercent,
			perGroupN, numGroups, betaScale, dsList);

			* build the name of the contrast data set;
		contrastDataSet = "contrasts" + "_" + simprefix;
		empiricalPowerDataSet = "empiricalPower" + "_" + simprefix; 

		* cleanup the SAS environment before we start;
		submit contrastDataSet;
		proc datasets;
			delete &contrastDataSet;
		run;
		endsubmit;

		* set the proc mixed call;
		if numGroups = 2 then mixedModMacroName="longit2Group5Rm";
		else mixedModMacroName="longit4Group5Rm";

		* set the covariance option for prox mixed;
		if covarCS then covariance='CS';
		else do;
			if covarCSH then covariance = 'CSH';
			else covariance = 'AR(1)';
		end;

		do i = 1 to nrow(dsList);
		* get the data set name;
		sdCurrent = dsList[i];
		*print sdCurrent;
		* call proc mixed;
		submit sdCurrent mixedModMacroName contrastDataSet covariance;
			ods exclude all;
			ods noresults;
			* fit the model and output p-value for test;
			ods output Contrasts=tmpSimContrasts;
			%&mixedModMacroName(&sdCurrent, &covariance);
			* append the results to the running data set;
			proc append base = &contrastDataSet data = tmpSimContrasts force; run;
			ods results;
			ods exclude none;
		endsubmit;
		end;

		* put the SAS environment back the way we found it;
		submit contrastDataSet empiricalPowerDataSet;
		data &contrastDataSet;
		  set &contrastDataSet;
		  reject = (probf < 0.05);
		run;
		proc freq data=&contrastDataSet;
			tables reject / out=&empiricalPowerDataSet;
		run;
		* stupid BS since the dynamic name only works on 64bit SAS;
		data temp;
			set &empiricalPowerDataSet;
			where reject = 1;
			keep percent;
		run;
		endsubmit;

		use temp;
		read all into empiricalPower;
		close temp;

		* remove the temp data set;
		call delete("work", "temp"); 

		* convert to a decimal from a percent;
		empiricalPower = empiricalPower / 100;
	finish;
								
	* load the parameter data into a matrix;
	paramNames={"targetPower" "monotone" "covarCS" "covarCSH" "covarAR1" 
				"missingPercent" "perGroupN" 
				"numGroups" "maxObservations" "betaScale"};
	use longitudinalParams;
  	read all var{"targetPower" "monotone" "covarCS" "covarCSH" "covarAR1" 
				"missingPercent" "perGroupN" 
				"numGroups" "maxObservations" "betaScale"} into paramList[colname=paramNames];

	*print paramList;
	
	call randseed(1546); 
	start = 1;
	resultNames = paramNames || "empiricalPower" || "empiricalElapsedTime";
	/*
	* Calculate empirical power for each 
	*/
	do i=1 to NROW(paramList);
		print ("Case: " + strip(char(i)));
		startTime = time();
		run longitEmpiricalPower(10000, 250, "outData", "simData",
							paramList[i,"monotone"],
							paramList[i,"covarCS"],
							paramList[i,"covarCSH"],
							paramList[i,"covarAR1"],
							paramList[i,"maxObservations"], 
							paramList[i,"missingPercent"],
							paramList[i,"perGroupN"], 
							paramList[i,"numGroups"],
							paramList[i,"betaScale"],
							empiricalPower);
		empiricalTime = time() - startTime;
		empiricalPowerResults = empiricalPowerResults // (empiricalPower || empiricalTime);

		if mod(i,10) = 0 | i = NROW(paramList) then do;
					* create final result set;
			results = paramList[start:i,] || empiricalPowerResults;
			* write power results to a data set;
			create tmpEmpirical from 
				results[colname=resultNames];
			append from results;
			close tmpEmpirical;
			free empiricalPowerResults;

		    dataSetName = "outData.empiricalPowerIter" + strip(char(start)) + "to" + strip(char(i));
			dsList = dsList // dataSetName;
			  /*
			  * Change the name of the data set using SAS
			  */
			  submit dataSetName;
				data &dataSetName;
					set tmpEmpirical;
				run;
			  endsubmit;
			start=i+1;
		end;
	end;
	

quit;

* append all of the data sets together;
data longitudinalEmpiricalPower;
	set 
	outData.empiricalPowerIter1to10
	outData.empiricalPowerIter11to20
	outData.empiricalPowerIter21to30
	outData.empiricalPowerIter31to40
	outData.empiricalPowerIter41to50
	outData.empiricalPowerIter51to60
	outData.empiricalPowerIter61to70
	outData.empiricalPowerIter71to80
	outData.empiricalPowerIter81to90
	outData.empiricalPowerIter91to100
	outData.empiricalPowerIter101to110
	outData.empiricalPowerIter111to120
	outData.empiricalPowerIter121to130
	outData.empiricalPowerIter131to140
	outData.empiricalPowerIter141to150
	outData.empiricalPowerIter151to160
	outData.empiricalPowerIter161to170
	outData.empiricalPowerIter171to180
	outData.empiricalPowerIter181to190
	outData.empiricalPowerIter191to200
	outData.empiricalPowerIter201to210
	outData.empiricalPowerIter211to220
	outData.empiricalPowerIter221to230
	outData.empiricalPowerIter231to240
	;
run;

* write the temporary empirical power data set to disk as a csv;
proc export data=longitudinalEmpiricalPower
   outfile="&OUT_DATA_DIR\longitudinalEmpiricalPower.csv"
   dbms=csv
   replace;
run;
