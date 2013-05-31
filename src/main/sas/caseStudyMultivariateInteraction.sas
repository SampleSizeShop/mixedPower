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


PROC IML SYMSIZE=1000 WORKSIZE=2000;
%INCLUDE "&POWERLIB_IML_FILE"/NOSOURCE2;
%INCLUDE "XMLUTILITIES.IML"/NOSOURCE2;

OPT_ON = {ORTHU UN HF GG BOX HLT PBT WLK MMETHOD UMETHOD MMETHOD};
* Specifying the option ORTHU in OPT_ON allows the program to provide;
* an orthonormal U matrix if one is not given by the user;
* This is the case for the following code;

* Create contrast matrices C and U (non-orthonormal) *;
P = 3;
Q = 4;
C = J(Q-1,1,1)||(-I(Q-1));
U = ( J(P-1,1,1)||(-I(P-1)) )`;

ALPHA = .01;

VARIANCE = 1;  
RHO = 0.4;
SIGMA = VARIANCE#(I(P)#(1-RHO) + J(P,P,RHO)); *Compound symmetry;
SIGSCAL = {1, 2};

ESSENCEX = I(Q); 
REPN = {5,10};
BETA = J(Q,P,0);
BETA[1,1] = 1;
BETASCAL = DO(0, 2.0 , 0.50);

/* MMETHOD = {4,4,4}; * Two moment null approximations + OBrien and Shieh (1992) 
					 noncentrality multiplier ON;

UCDF = {4,2,2,4};  * UN and Box (4):
					 Exact via Davies' algorithm (1980), as in Muller, 
					 Edwards, Simpson, and Taylor (2007). If exact fails, 
					 then switch to approximation 2, MEST (2007);
                   * HF and GG (2): 
					 Muller, Edwards, Simpson, and Taylor (2007) 
					 approximation;
*/

* Output full precision ;
ROUND = 15;

RUN POWER;

/* write the data to an XML file */
TEST_LIST = {'unirep' 'unirepBox' 'unirepGG' 'unirepHF' 'wl' 'pbt' 'hlt'};
filename out "&OUTPUT_DATA_DIRECTORY\TestConditionalMultivariateInteraction.xml";
RUN powerResultsToXML(out, _HOLDPOWER, TEST_LIST, 0);

QUIT;




* define the mixed model fitting macro; 
* this must contain a 'by setID' statement, but can otherwise;
* be defined as needed by the model;
%macro fitMixedModelShort(datasetName);
	proc mixed data=&datasetName;
		model y = A1 A2 A3 B1 B2 B3 C1 C2 C3 D1 D2 D3 / noint solution ddfm=KR;
		random int / subject=patID;
		by setID;
		contrast "Time x Treatment Interaction" 
		A1 1 A2 -1 A3 0 B1 -1 B2 1 B3 0 C1 0 C2 0 C3 0 D1 0 D2 0 D3 0,
		A1 1 A2 0 A3 -1 B1 -1 B2 0 B3 1 C1 0 C2 0 C3 0 D1 0 D2 0 D3 0,
		A1 1 A2 -1 A3 0 B1 0 B2 0 B3 0 C1 -1 C2 1 C3 0 D1 0 D2 0 D3 0,
		A1 1 A2 0 A3 -1 B1 0 B2 0 B3 0 C1 -1 C2 0 C3 -1 D1 0 D2 0 D3 0,
		A1 1 A2 -1 A3 0 B1 0 B2 0 B3 0 C1 0 C2 0 C3 0 D1 -1 D2 1 D3 0,
		A1 1 A2 0 A3 -1 B1 0 B2 0 B3 0 C1 0 C2 0 C3 0 D1 -1 D2 0 D3 1
		;
	run;
%mend;

/*
* Generate the data sets
*/
proc iml;
	%INCLUDE "&MODULES_DIR\simulateMixedModel.sxs"/NOSOURCE2;
	%INCLUDE "&MODULES_DIR\calculatePowerKenwardRoger.sxs"/NOSOURCE2;

* Create contrast matrices C and U (non-orthonormal) *;
P = 3;
Q = 4;
C = J(Q-1,1,1)||(-I(Q-1));
U = ( J(P-1,1,1)||(-I(P-1)) )`;
* convert to a mixed model contrast;
CFull = C@U`;
 
* design matrix for means;
ESSENCEX = I(P*Q);
REPN = 10;
X = J(REPN,1,1)@ESSENCEX;
ISU = DO(1,Q*REPN,1)@J(1,3,1);
*print ISU;

X = ISU` || X ;
*print X;
XFullColNames={"PatId" "A1" "A2" "A3" "B1" "B2" "B3" "C1" "C2" "C3" "D1" "D2" "D3"}; 
XModelColNames = {"A1" "A2" "A3" "B1" "B2" "B3" "C1" "C2" "C3" "D1" "D2" "D3"};
mattrib X colname=XFullColNames; 
Xmodel = X[,XModelColNames];

* beta (means) ;
BETA = J(Q*P,1,0);
BETA[1,1] = 1;
BETASCAL = {1}; *DO(0, 2.0 , 0.50);

* covariance matrix - compound symmetric within repeated measures;
VARIANCE = 1;  
RHO = 0.4;
SIGMA = VARIANCE#(I(P)#(1-RHO) + J(P,P,RHO)); *Compound symmetry;
print SIGMA;
SIGMAS = I(Q*REPN)@SIGMA;

thetaNull = J(6,1,0);
* Type I error rate;
ALPHA = .05;

simlib= "outData";
simprefix = "multInt";

*call calculateEmpiricalPowerConditional(10000, 1000,  
  simlib, simprefix, "fitMixedModelShort",
  X, XFullColNames, XModelColNames, Beta, SigmaS,
  powerResults);

* now call the approximation;
power = calculatePowerKenwardRoger(Xmodel, Beta, CFull, SigmaS, thetaNull, alpha, Q*REPN);
print power;
print powerResults;
quit;





