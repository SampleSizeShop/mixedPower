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
		model y = A B C D / noint ddfm=KR;
		by setID;
		contrast "trt" A 1 B -1 C 0 D 0,
					A 1 B 0 C -1 D 0,
					A 1 B 0 C 0 D -1;
	run;
%mend;

%macro fitGLM(datasetName);
	proc glm data=&datasetName noprint;
		model y = A B C D / noint;
		by setID;
		contrast "trt" A 1 B -1 C 0 D 0,
					A 1 B 0 C -1 D 0,
					A 1 B 0 C 0 D -1;
	run;
%mend;

/*
* Generate the data sets
*/
proc iml;
	%INCLUDE "simulateMixedModel.sxs"/NOSOURCE2;
	%INCLUDE "calculatePowerKenwardRoger.sxs"/NOSOURCE2;
Xessence = I(4);

X = Xessence@J(5,1,1);
XFullColNames={"A" "B" "C" "D"}; 
XModelColNames = XFullColNames;
mattrib X colname=XFullColNames; 
Xmodel = X[,XModelColNames];

Beta = {1,0,0,0};

SigmaS = 2#I(20);
C = {1 -1 0 0, 1 0 -1 0, 1 0 0 -1};
thetaNull = {0,0,0};
alpha = {0.05};
simlib= "outData";
simprefix = "simCaseStudy03";

call calculateEmpiricalPowerConditional(10000, 1000,  
  simlib, simprefix, "fitMixedModelShort",
  X, XFullColNames, XModelColNames, Beta, SigmaS,
  empiricalPower);

* now call the approximation;
calculatedPower = calculatePowerKenwardRoger(Xmodel, Beta, C, SigmaS, thetaNull, alpha, nrow(X));
print calculatedPower[colname={"Fcrit" "Fobserved" "ndf" "ddf" "omegaP" "power"}];

powerTable = calculatedPower[,6] || empiricalPower;
print powerTable[colname={"Calculated Power" "Empirical Power"}];

quit;





