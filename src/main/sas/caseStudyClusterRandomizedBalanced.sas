/*
*
* Case Study 4: Balanced data with 40 outcomes
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
		random int / subject=clusterID;
		by setID;
		contrast "trt" A 1 B -1;
	run;
%mend;

/*
* Generate the data sets
*/
proc iml;
	%INCLUDE "simulateMixedModel.sxs"/NOSOURCE2;
	%INCLUDE "calculatePowerKenwardRoger.sxs"/NOSOURCE2;
Xessence = {
	1 0 ,
	0 1
};

X = Xessence@J(200,1,1);
ISU = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}@J(40,1,1);

X = ISU || X ;
XFullColNames={"clusterId" "A" "B"}; 
XModelColNames = {"A" "B" };
mattrib X colname=XFullColNames; 
Xmodel = X[,XModelColNames];

Beta = {1,0};

SigmaI = 2#(J(40,40,1)*0.3+I(40)*0.7);
SigmaS = I(10)@SigmaI;
print SigmaI;
C = {1 -1};
thetaNull = {0};
alpha = {0.05};
simlib= "outData";
simprefix = "simCaseStudy04";

call calculateEmpiricalPowerConditional(10000, 1000,  
  simlib, simprefix, "fitMixedModelShort",
  X, XFullColNames, XModelColNames, Beta, SigmaS,
  powerResults);

* now call the approximation;
power = calculatePowerKenwardRoger(Xmodel, Beta, C, SigmaS, thetaNull, alpha, 10);
print power;
print powerResults;
quit;





