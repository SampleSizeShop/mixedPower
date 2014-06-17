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

proc iml;
	foo = I(3);
	print foo;
quit;

/*
* Calculate empirical power for the 4 group, cluster randomized trials
*/
proc iml;
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
		return(matrix);
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
		return(matrix);
	finish;

	/* 
	* generate a heterogeneous compound symmetric covariance
	* sigmaList - an Nx1 vector of sigma values
	* rho - the correlation between observations
	*/
	start heterogeneousCSMatrix(sigmaList,rho, matrix) {
		size = nrow(sigmaList);
		matrix = J(size,size,0);
		do r = 1 to size;
			do c = 1 to size;
		      if (r==c) {
		        matrix[r,c] = sigmaList[r]*sigmaList[c]
		      } else {
		        matrix[r,c] = sigmaList[r]*sigmaList[c] * rho
		      }
			end;
		end;
		return(matrix);
	finish;


	/*
	* Generate data for a longitudinal design with the specified parameters
	*/
	start generateData(replicates, setSize, simLib, simPrefix,
		monotone, maxObservations, missingPercent,
		perGroupN, numGroups, betaScale);

		* determine the number of complete and incomplete sampling units;
		numComplete = floor(perGroupN * (1 - missingPercent);
		numIncomplete = perGroupN - numComplete;

		* build the within design matrix for complete cases and incomplete cases;
		designWithinComplete = I(maxObservations);
		X = J(numComplete,1,1) @ I(numGroups) @ designWithinComplete; 
		if numIncomplete != 0 then;
			if missingType = "monotone" then do;
			    designWithinIncomplete = I(maxObservations)[1:(maxObservations-2),];
			else;
				designWithinIncomplete = I(maxObservations)[{1,3,5},];
			end;
			X = X // J(numIncomplete,1,1) @ I(numGroups) @ designWithinIncomplete; 
		end;

		* build the beta matrix;
		Beta = {1} // J(((numGroups * maxObservations) - 1), 1, 0);

		* create the Sigma matrix for the complete cases and incomplete cases; 
		rho = 0.04;
		if covariance = "CS" then;
			heterogeneousCSMatrix(J(maxObservations,1,1), rho, SigmaComplete);
		else;
			if covariance = "CSH" then;
				heterogeneousCSMatrix({1.0.5,0.3,0.1,0.1}, rho, SigmaComplete);
			else;
				ar1Matrix(maxObservations, rho, SigmaComplete);
			end;
		end;
		if numIncomplete != 0 then;
			SigmaIncomplete = designWithinIncomplete * SigmaComplete * designWithinIncomplete`;
		end;

		* build the data set column names;
		if numGroups = 2 then do;
			XFullColNames={"subjectId" "trt1_rep1" "trt1_rep2" "trt1_rep3" "trt1_rep4" "trt1_rep5"
					"trt2_rep1" "trt2_rep2" "trt2_rep3" "trt2_rep4" "trt2_rep5"}; 
		else do;
			XFullColNames={"subjectId" "trt1_rep1" "trt1_rep2" "trt1_rep3" "trt1_rep4" "trt1_rep5"
			"trt2_rep1" "trt2_rep2" "trt2_rep3" "trt2_rep4" "trt2_rep5" 
			"trt3_rep1" "trt3_rep2" "trt3_rep3" "trt3_rep4" "trt3_rep5" 
			"trt4_rep1" "trt4_rep2" "trt4_rep3" "trt4_rep4" "trt4_rep5" }; 
		end;

		* determine the number of sets;
		numSets = replicates / setSize;

		* create a vector of means (all 0's) for the random errors;
		mu = J(maxObservations,1,0);

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
		    E = E // randnormal(numComplete * numGroups, mu, SigmaComplete)
			idList = idList // (T(isu:isu+numComplete-1) @ J(maxObservations,1,1));
			isu = isu + numComplete;

			if numIncomplete > 0 then;
			    EIncomplete = randnormal(numIncomplete * numGroups, mu, SigmaComplete)
				idList = idList // (T(isu:isu+numComplete-1) @ J(maxObservations,1,1));
			end;
			* form the responses;
			Y = X * Beta + E;
			* build the data set;
			block = J(nrow(Y),1,setId) || Y || E || idList || X;
			setId = setId + 1;

			* append the block to the data set;
			output = output // block;
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
		  *//*
		  submit dataSetName;
			data &dataSetName;
				set temp;
			run;
		  endsubmit;
		  */
		  * remove the temp data set;
		  call delete("work", "temp"); 
		end;

	finish;


	* wrapper function to generate appropriate X and Sigma;
	start longitEmpiricalPower(covariance, missingType, missingPercent, 
							perGroupN, numGroups, maxObservations, betaScale,
							empiricalPower);

		* generate data for the design;
		generateData(10000, 1000, "work", "simData",
			monotone, maxObservations, missingPercent,
			perGroupN, numGroups, betaScale);


	finish;
								
	* load the parameter data into a matrix;
	paramNames={"targetPower" "covariance" "missingType" "missingPercent" "perGroupN" 
				"numGroups" "maxObservations" "betaScale"};
	use longitudinalParams;
  	read all into paramList[colname=paramNames];

	do i=1 to NROW(paramList);
		print ("Case: " + strip(char(i)));

		call longitEmpiricalPower(paramList[i,"covariance"],
							paramList[i,"missingType"],
							paramList[i,"missingPercent"],
							paramList[i,"perGroupN"], 
							paramList[i,"numGroups"],
							paramList[i,"maxObservations"], 
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
