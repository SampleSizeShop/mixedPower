/*
*
* Case Study 1: Two sample t-test.
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
		by setID;
		contrast "treatment effect" A 1 B -1;
	run;
%mend;

/*
* Generate the data sets
*/
proc iml;
	%INCLUDE "simulateMixedModel.sxs"/NOSOURCE2;
	%INCLUDE "calculatePowerKenwardRoger.sxs"/NOSOURCE2;
	%INCLUDE "C:\KeithMullerSoftware\power\Iml\POWERLIB21.IML" /NOSOURCE2;
essenceX = {
	1 0 ,
	0 1
};

X = essenceX@J(10,1,1);

XFullColNames={"A" "B"}; 
XModelColNames = {"A" "B"};

Beta = {1,0};

SigmaS = 2#I(20);
C = {1 -1};
thetaNull = {0};
alpha = {0.05};
simlib= "outData";
simprefix = "simCaseStudy01";

*call calculateEmpiricalPowerConditional(10000, 1000,  
  simlib, simprefix, "fitMixedModelShort",
  X, XFullColNames, XModelColNames, Beta, SigmaS,
  powerResults);

* now call the approximation;
power = calculatePowerKenwardRoger(X, Beta, C, SigmaS, thetaNull, alpha, nrow(X));
print power;
print powerResults;

sigma = {2};
REPN = {10};
OPT_ON = {HLT };
OPT_OFF= {COLLAPSE C GG HF UN BOX};
ROUND = 15;
RUN POWER;

quit;





