/*
*
* Case Study 5: Unalanced data with 5:10:10:5:5 outcomes
* per independent sampling unit.
*
* Fixed imbalance
*
*
*
*/

%include "common.sas";

* define the mixed model fitting macro; 
* this must contain a 'by setID' statement, but can otherwise;
* be defined as needed by the model;
%macro fitMixedModelShort(datasetName);
	proc mixed data=&datasetName;
		model y = A B / noint solution ddfm=KR;
		random int / subject=clusterID;
		by setID;
		contrast "trt" A 1 B -1;
	run;
%mend;

/*
* Generate the data sets
*/
proc iml;
	%INCLUDE "&MODULES_DIR\simulateMixedModel.sxs"/NOSOURCE2;
	%INCLUDE "&MODULES_DIR\calculatePowerKenwardRoger.sxs"/NOSOURCE2;
Xessence = {
	1 0 ,
	0 1
};

X = Xessence@J(35,1,1);
ISU = J(5,1,1) // J(10,1,2) // J(10,1,3) // J(5,1,4) // J(5,1,5) //
		J(5,1,6) // J(10,1,7) // J(10,1,8) // J(5,1,9) // J(5,1,10);

X = ISU || X ;
XFullColNames={"clusterId" "A" "B"}; 
XModelColNames = {"A" "B" };
mattrib X colname=XFullColNames; 
Xmodel = X[,XModelColNames];

Beta = {1,0};

SigmaSmall = 2#(J(5,5,1)*0.3+I(5)*0.7);
SigmaBig = 2#(J(10,10,1)*0.3+I(10)*0.7);
SigmaS = block(SigmaSmall,SigmaBig,SigmaBig,SigmaSmall,SigmaSmall,
				SigmaSmall,SigmaBig,SigmaBig,SigmaSmall,SigmaSmall);
print SigmaS;
C = {1 -1};
thetaNull = {0};
alpha = {0.05};
simlib= "outData";
simprefix = "simCaseStudy05";

call calculateEmpiricalPowerConditional(10000, 1000,  
  simlib, simprefix, "fitMixedModelShort",
  X, XFullColNames, XModelColNames, Beta, SigmaS,
  powerResults);

* now call the approximation;
power = calculatePowerKenwardRoger(Xmodel, Beta, C, SigmaS, thetaNull, alpha, 10);
print power;
print powerResults;
quit;





