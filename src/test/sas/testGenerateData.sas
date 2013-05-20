
libname outData "../output/datasets";

proc iml;
	%INCLUDE "simulateMixedModel.sxs"/NOSOURCE2;
	
start myfun(A, cols, sub);
*mattrib A colname=cols; 
print A;
foo1 = A[,sub];
print foo1;
finish;

Xessence = {
	1 0 0 0,
	0 1 0 0,
	0 0 1 0,
	0 0 0 1
};

ISU = {1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10};

X = Xessence@J(5,1,1);

X = ISU || X ;
XFullColNames={"PatId" "A" "B" "C" "D"}; 
XModelColNames = {"A" "B" "C" "D"};

Beta = {1,0,0,0};

SigmaI = 2#{
	1 0.2 ,
	0.2 1 
};
SigmaS = I(10)@SigmaI;

datasetPrefix= "outData.mysim";

call genDataFullCovarianceNoMissing(10000, 1000, datasetPrefix, 
	X, XFullColNames, XModelColNames, Beta, SigmaS);


quit;
