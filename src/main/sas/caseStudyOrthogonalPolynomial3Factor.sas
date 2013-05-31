/*
*
* Case Study: 3 factor design with orthogonal polynomial contrasts.
*
* The study design for Example 1 is a balanced, two-group design. 
* We calculate power for a two-sample t-test comparing the mean 
* responses between the two groups. The example is based on the results in
* Muller, K. E., & Benignus, V. A. (1992). 
* Neurotoxicology and teratology, 14(3), 211-219.
*
*/

%include "common.sas";

PROC IML SYMSIZE=1000 WORKSIZE=2000;
%INCLUDE "&POWERLIB_IML_FILE"/NOSOURCE2;
%INCLUDE "XMLUTILITIES.IML"/NOSOURCE2;
RESET FUZZ NOAUTONAME FW=6 LINESIZE=80;

ALPHA = .05;

* Choose dimensions of design *;
GA = 3; * =# groups for between factor A *;
GB = 3; * =# groups for between factor B *;
GC = 3; * =# groups for between factor C *;
TD = 3; * =#Times for within factor D *;
TE = 3; * =#Times for within factor E *;
TF = 3; * =#Times for within factor F *;

P = TD#TE#TF;
Q = GA#GB#GC;
ESSENCEX = I(Q);
BETA = J(Q,P,0);
BETA[1,1] = 1;
REPN = DO(2,12,2);
SIGMA = DIAG(DO(1,P,1)); * Variances are 1,2,3,...p *;

* Get orthonormal U matrices *;
CALL UPOLY3 ( (1:TD),"D", (1:TE),"E",  (1:TF),"F",
		          UD,UDLBL,   UE,UELBL,    UF,UFLBL, 
                 UDE,UDELBL, UDF,UDFLBL,  UEF,UEFLBL,  UDEF,UDEFLBL );
 
* Get orthonormal C matrices *;
CALL UPOLY3 ((1:GA),"A" , (1:GB),"B" , (1:GC),"C",
		         U1,CALBL,    U2,CBLBL,    U3,CCLBL,
                U12,CABLBL, U13,CACLBL,   U23,CBCLBL,  U123,CABCLBL);

BETASCAL = {9 18 27};
ROUND = 15;

OPT_ON = {UN HF GG BOX  HLT PBT WLK MMETHOD  UMETHOD MMETHOD NOPRINT};
OPT_OFF = {COLLAPSE};
BUG=" ";

C = U1`;
U = UD;
RUN POWER;
HOLDA=HOLDA//_HOLDPOWER;
PRINT / "AxD";
PRINT HOLDA[COLNAME=_HOLDPOWERLBL ROWNAME=BUG];

C = U12`;
U = UDE;
RUN POWER;
HOLDABDE = HOLDABDE//_HOLDPOWER;
PRINT / "AxB x DxE Interaction";
PRINT HOLDABDE[COLNAME=_HOLDPOWERLBL ROWNAME=BUG];

C = U123`;
U = UDEF;
RUN POWER;
HABCDEF = HABCDEF//_HOLDPOWER;
PRINT / "AxBxC x DxExF Interaction";
PRINT HABCDEF[COLNAME=_HOLDPOWERLBL ROWNAME=BUG];

HOLDALL = HOLDA//HOLDABDE//HABCDEF;
PRINT (NROW(HOLDALL));
PRINT HOLDALL[COLNAME=_HOLDPOWERLBL];
/* write the data to an XML file */
TEST_LIST = {'unirep' 'unirepBox' 'unirepGG' 'unirepHF' 'wl' 'pbt' 'hlt'};
filename out "&OUTPUT_DATA_DIRECTORY\TestConditionalOrthogonalPolynomial3Factor.xml";
RUN powerResultsToXML(out, HOLDALL, TEST_LIST, 0);

QUIT;


* define the mixed model fitting macro; 
* this must contain a 'by setID' statement, but can otherwise;
* be defined as needed by the model;
%macro fitMixedModel(datasetName);
	proc mixed data=&datasetName;
		model y = A B / noint solution ddfm=KR;
		by setID;
		contrast "treatment effect" A 1 B -1;
	run;
%mend;

/*
* We perform the following steps:
* 1. Calculate empirical power for a mixed model with Kenward-Roger ddfm
* 2. Calculate power using our noncentral F approximation
* 3. Calculate empirical power for a general linear model
* 4. Calculate power using power methods described by Muller et al.
*/
proc iml;
	%INCLUDE "&MODULES_DIR\simulateMixedModel.sxs"/NOSOURCE2;
	%INCLUDE "&MODULES_DIR\calculatePowerKenwardRoger.sxs"/NOSOURCE2;
	%INCLUDE "&MODULES_DIR\POWERLIB21.IML" /NOSOURCE2;

* design matrix for means;
essenceX = {
	1 0 ,
	0 1
};
X = essenceX@J(10,1,1);

* Column names;
XFullColNames={"A" "B"};
* Names of columns involved in the model; 
XModelColNames = {"A" "B"};

* mean difference between the groups;
BaseBeta = {1,0};

* covariance matrix.  assumes independence of all observations;
BaseSigmaS = I(20);

* between participant contrast, testing difference between the groups;
C = {1 -1};

* null hypothesis;
thetaNull = {0};

* Type I error rate;
alpha = {0.05};

* prefixes for simulated data sets;
simlib= "outData";
simprefix = "caseStudy01";

* list defining the mean differences between the groups;
betascale = {0 1};
* do(0,2,1);
* list defining the covariance values;
sigmascale = {0.32 1.00 2.05};

* calculate power for several mean differences and covariance values;
do bIndex = 1 to ncol(betascale);
  beta = betaScale[bIndex]#basebeta;
  do sIndex = 1 to ncol(sigmaScale);
    sigmaS = sigmaScale[sIndex]#baseSigmaS;

    * get simulated power for the mixed model;
    call calculateEmpiricalPowerConditional(10000, 1000,  
            simlib, simprefix, "fitMixedModel",
            X, XFullColNames, XModelColNames, Beta, SigmaS,
            powerResults);

    * now call the approximation;
    power = calculatePowerKenwardRoger(X, Beta, C, SigmaS, thetaNull, alpha, nrow(X));

	results = results // (betaScale[bindex] || sigmaScale[sIndex] || power[,6] || powerResults);
  end;
end;

print results[colname={"Beta-scale" "Sigma-scale" "Calculated Power" "Empirical Power"}];

/*
sigma = {2};
REPN = {10};
OPT_ON = {HLT };
OPT_OFF= {COLLAPSE C GG HF UN BOX};
ROUND = 15;
RUN POWER;
*/
quit;





