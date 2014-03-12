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
		model y = trt1 trt2 trt3 trt4 / noint solution ddfm=kr;
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
		model y = trt1 trt2 / noint solution ddfm=kr;
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
	set clusterRandomizedParams;
	where clusterSize=50 and perGroupN=40;
run;

data clusterRandomizedParams;
	set clusterRandomizedParams(firstobs=10 obs=18);
run;

ods graphics off;
/*
* Calculate empirical power for the 4 group, cluster randomized trials
*/
proc iml;

	start genData(replicates, setSize, simlibname, simPrefix, 
				numGroups, perGroupN, clusterSize, 
				deletePercent, sigmaSq, icc, beta, 
				XColNames, dsList);

		datasetPrefix = simlibname + "." + simprefix;

		* determine the number of sets;
		numSets = replicates / setSize;

		* calculate the number of complete cases;
		numComplete = perGroupN / 2;
		numIncomplete = perGroupN - numComplete;
		* calculate the size of complete and incomplete cases;
		incompleteSize = floor(clusterSize*(1-deletePercent));

		* get sigmas for complete and incomplete cases;
		SigmaComplete = sigmaSq#(J(clusterSize,clusterSize,1)*icc+I(clusterSize)*(1-icc));
		SigmaIncomplete = sigmaSq#(J(incompleteSize,incompleteSize,1)*icc+I(incompleteSize)*(1-icc));
	 
		FullColNames = {"setID" "Y" "E"} || XColnames;
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
		do setNum = 1 to numSets;
		  do repNum = 1 to setSize;

		  	free X;
			free E;
			free Y;
			* generate data for a single data set;
		  	isu = 1;
			do grp = 1 to numGroups;
				
				XComplete = J(clusterSize,1,1) @ I(numGroups)[grp,];
				muComplete = J(clusterSize,1,0);
				XIncomplete = J(incompleteSize,1,1) @ I(numGroups)[grp,];
				muIncomplete = J(incompleteSize,1,0);
				Xpair = XComplete // XIncomplete;
				* assumes equal numbers of complete and incomplete;
				do pair=1 to numComplete;
					Epair = randnormal(1, muComplete, SigmaComplete) || 
						randnormal(1, muIncomplete, SigmaIncomplete);
					* transpose into long format;
					ETpair = Epair`;
					* calculate Y;
					Ypair = Xpair * Beta + ETpair;

					X = X // ((J(clusterSize,1,isu) // J(incompleteSize,1,isu+1)) || XPair);
					E = E // ETpair;
					Y = Y // Ypair;
					isu = isu + 2;
				end;
			end;


			* build the data set;
			block = (J(NROW(Y),1,repNum + ((setNum-1)*setSize)) || Y || E || X);
			mattrib block colname=(FullColNames); 

			* append the full X matrix to the data set;
			output = output // block;
		  end;
		  * write to disk;
		  startIterNum = (setNum-1)*setSize+1;
		  startIter = char(startIterNum);
		  endIterNum = startIterNum + setSize - 1;
		  endIter = char(endIterNum); 
		  dataSetName = datasetPrefix + "Iter" + strip(startIter) + 
				"to" + strip(endIter);

		  dsList = dsList // dataSetname;

		  /*
		  * Write a temporary data set since I
		  * can't convince IML to use a dynamic name here
		  */
		  create temp from output[colname=FullColNames];
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
	start clusterEmpiricalPower(numGroups, perGroupN, clusterSize, 
								deletePercent, sigmaSq, icc, betaScale,
								empiricalPower);

		if numGroups = 2 then do;
			XFullColNames={"clusterId" "trt1" "trt2"}; 
			* define beta;
			Beta = betaScale*{1,0};
			macroName = "clustered2Group";
		end;
		else do;
			XFullColNames={"clusterId" "trt1" "trt2" "trt3" "trt4"}; 
			* define beta;
			Beta = betaScale*{1,0,0,0};
			macroName = "clustered4Group";
		end;

		simlib= "outData";
		simprefix = "clusterLargeEx";
		replicates=10000;
		blockSize = 100;

		call genData(replicates, blockSize, simlib, simPrefix, 
				numGroups, perGroupN, clusterSize, 
				deletePercent, sigmaSq, icc, beta, 
				XFullColNames, dsList);

		* build the name of the contrast data set;
		contrastDataSet = "contrasts" + "_" + simprefix; 
		empiricalPowerDataSet = "empiricalPower" + "_" + simprefix; 

		* cleanup the SAS environment before we start;
		submit contrastDataSet;
			proc datasets;
				delete &contrastDataSet;
			run;
  		endsubmit;

		do i = 1 to nrow(dsList);
		    * get the data set name;
			sdCurrent = dsList[i];
		    print sdCurrent;
			* call proc mixed;
		    submit sdCurrent macroName contrastDataSet;
				ods exclude all;
				ods noresults;
				* fit the model and output p-value for test;
				ods output Contrasts=tmpSimContrasts;
				%&macroName(&sdCurrent);
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
							paramList[i,"betaScale"],
							empiricalPower);
							print(empiricalPower);
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
proc export data=clusterRandomizedEmpirical
   outfile="&OUT_DATA_DIR\clusterLargeEmpirical.csv"
   dbms=csv
   replace;
run;
