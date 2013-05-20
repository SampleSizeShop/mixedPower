/*
*
* Case Study: Two sample t-test.
*
* The study design for Example 1 is a balanced, two-group design. 
* We calculate power for a two-sample t-test comparing the mean 
* responses between the two groups. The example is based on the results in
* Muller, K. E., & Benignus, V. A. (1992). 
* Neurotoxicology and teratology, 14(3), 211-219.
*
*/

%include "common.sas";

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





