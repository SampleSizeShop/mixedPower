/*
* Generated SAS code for simulating a linear mixed model
*
* Author: Sarah Kreidler
* Date: 8/12/2013
*/
{r setup, echo=FALSE}
matrixToIML <- function(name, m) {
  cat(name,"= {\n")
  rowCount = nrow(m)
  for(i in 1:rowCount) {
    cat("\t",m[i,])
    if (i < rowCount) {
      cat(",\n")
    }
  }
  cat("\n}\n")
}

PROC IML SYMSIZE=1000 WORKSIZE=2000;

%INCLUDE "calculatePowerKenwardRoger.sxs"/NOSOURCE2;

Xessence = {

```
X = {
	 1 0 0 0,
	 0 1 0 0,
	 0 0 1 0,
	 0 0 0 1
}
```

};

X = Xessence@J(10,1,1);

C = {1 -1};

SigmaS = I(NROW(X));

Beta = {
	1,
	0
};

thetaNull = {0};
alpha = {0.05};

print X;
print sigmaS;
print Beta;
print C;

do i = 0 to 3;
	power = power // calculatePowerKenwardRoger(X, i*Beta, C, SigmaS, thetaNull, alpha);
end;
print power;

quit;


/*
*
* POwerlib equivalent
*
*/
PROC IML SYMSIZE=1000 WORKSIZE=2000;
%INCLUDE "C:\KeithMullerSoftware\power\Iml\POWERLIB21.IML" /NOSOURCE2;

* Define inputs to power program;
ALPHA = 0.05;
SIGMA = {1};
SIGSCAL = {1};

ESSENCEX = I(2);
REPN = { 10 };

BETA = {1 0}`;
BETASCAL = DO(0,3,1);
C = {1 -1};

OPT_OFF= {C U};
ROUND = 15;
RUN POWER;

QUIT;
