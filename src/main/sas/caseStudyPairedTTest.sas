/*
*
* Case Study 1: Balanced data with 2 outcomes
* per independent sampling unit.
*
*
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
		random int / subject=patID;
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

X = Xessence@J(10,1,1);
ISU = {1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10};
X = ISU || X ;
XFullColNames={"PatId" "A" "B"}; 
XModelColNames = {"A" "B" };
mattrib X colname=XFullColNames; 
Xmodel = X[,XModelColNames];

Beta = {1,0};

SigmaI = 2#{
	1 0.3 ,
	0.3 1 
};
SigmaS = I(10)@SigmaI;
C = {1 -1};
thetaNull = {0};
alpha = {0.05};
simlib= "outData";
simprefix = "simCaseStudy02";

call calculateEmpiricalPowerConditional(10000, 1000,  
  simlib, simprefix, "fitMixedModelShort",
  X, XFullColNames, XModelColNames, Beta, SigmaS,
  powerResults);

* now call the approximation;
power = calculatePowerKenwardRoger(Xmodel, Beta, C, SigmaS, thetaNull, alpha, 10);
print power;
print powerResults;
quit;





