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
		model y = A B / noint solution ddfm=kr;
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

X = Xessence@J(50,1,1);

*ISU = J(10,1,1) // J(10,1,2) // J(10,1,3) // J(5,1,4) // J(5,1,5) //
		J(10,1,6) // J(10,1,7) // J(10,1,8) // J(5,1,9) // J(5,1,10);
ISU = J(5,1,1) // J(5,1,2) // J(5,1,3) // J(5,1,4) // J(5,1,5)
	// J(5,1,6) // J(5,1,7) // J(5,1,8) // J(5,1,9) // J(5,1,10)
	// J(5,1,11) // J(5,1,12) // J(5,1,13) // J(5,1,14) // J(5,1,15)
	// J(5,1,16) // J(5,1,17) // J(5,1,18) // J(5,1,19) // J(5,1,20);
X = ISU || X ;
XFullColNames={"clusterId" "A" "B"}; 
XModelColNames = {"A" "B" };
mattrib X colname=XFullColNames; 
Xmodel = X[,XModelColNames];

*print X;
Beta = {1,0};

SigmaSmall = 2#(J(5,5,1)*0.3+I(5)*0.7);
SigmaBig = 2#(J(10,10,1)*0.3+I(10)*0.7);
SigmaS = I(20) @ SigmaSmall;
*print SigmaS;
C = {1 -1};
thetaNull = {0};
alpha = {0.05};
simlib= "outData";
simprefix = "balanced";

call calculateEmpiricalPowerConditional(10000, 1000,  
  simlib, simprefix, "fitMixedModelShort",
  X, XFullColNames, XModelColNames, Beta, SigmaS,
  empiricalPower);

* now call the approximation;
power = calculatePowerKenwardRoger(Xmodel, Beta, C, SigmaS, thetaNull, alpha, 10);
print power;
print empiricalPower;
quit;


******** A > 1******;
* define the mixed model fitting macro; 
* this must contain a 'by setID' statement, but can otherwise;
* be defined as needed by the model;
%macro fitMixedA2(datasetName);
	proc mixed data=&datasetName;
		model y = A B C / noint solution ddfm=kr;
		random int / subject=clusterID;
		by setID;
		contrast "trt ABC" A 1 B -1, A 1 C -1;
	run;
%mend;

/*
* Generate the data sets
*/
proc iml;
	%INCLUDE "&MODULES_DIR\simulateMixedModel.sxs"/NOSOURCE2;
	%INCLUDE "&MODULES_DIR\calculatePowerKenwardRoger.sxs"/NOSOURCE2;
Xessence = I(3);

X = Xessence@J(50,1,1);

*ISU = J(10,1,1) // J(10,1,2) // J(10,1,3) // J(5,1,4) // J(5,1,5) //
		J(10,1,6) // J(10,1,7) // J(10,1,8) // J(5,1,9) // J(5,1,10);
ISU = J(5,1,1) // J(5,1,2) // J(5,1,3) // J(5,1,4) // J(5,1,5)
	// J(5,1,6) // J(5,1,7) // J(5,1,8) // J(5,1,9) // J(5,1,10)
	// J(5,1,11) // J(5,1,12) // J(5,1,13) // J(5,1,14) // J(5,1,15)
	// J(5,1,16) // J(5,1,17) // J(5,1,18) // J(5,1,19) // J(5,1,20)
	// J(5,1,21) // J(5,1,22) // J(5,1,23) // J(5,1,24) // J(5,1,25)
	// J(5,1,26) // J(5,1,27) // J(5,1,28) // J(5,1,29) // J(5,1,30);
X = ISU || X ;
XFullColNames={"clusterId" "A" "B" "C"}; 
XModelColNames = {"A" "B" "C"};
mattrib X colname=XFullColNames; 
Xmodel = X[,XModelColNames];

*print X;
Beta = {1,0,0};

SigmaSmall = 2#(J(5,5,1)*0.3+I(5)*0.7);
SigmaBig = 2#(J(10,10,1)*0.3+I(10)*0.7);
SigmaS = I(30) @ SigmaSmall;
*print SigmaS;
C = {1 -1 0, 1 0 -1};
thetaNull = {0};
alpha = {0.05};
simlib= "outData";
simprefix = "balanced3Grp";

call calculateEmpiricalPowerConditional(10000, 1000,  
  simlib, simprefix, "fitMixedA2",
  X, XFullColNames, XModelColNames, Beta, SigmaS,
  empiricalPower);

* now call the approximation;
power = calculatePowerKenwardRoger(Xmodel, Beta, C, SigmaS, thetaNull, alpha, 10);
print power;
print empiricalPower;
quit;



/*
data test;
set outdata.balanced3grpiter1to1000;
where setID = 1;
run;

proc mixed data=test;
	model y = A B C / noint solution ddfm=kr;
	random int / subject=clusterID;
	contrast "trt ABC" A 1 B -1, A 1 B -1;
run;*/
