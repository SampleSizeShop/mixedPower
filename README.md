The mixedPower Package
=========================

The mixedPower package for R (>3.0.0) calculates power for the 
linear mixed model for longitudinal data. The power calculations
may be used for designs with or without anticipated missing data.  

The package provides companion code for the manuscript:

Kreidler, S. M., Muller, K. E., & Glueck, D. H. 
A Power Approximation for Longitudinal Studies Using the 
Kenward and Roger Wald Test in the Linear Mixed Model, In review.

### Instructions for replicating the manuscript results 

The results in the above manuscript were produced using R version 3.0.2. To reproduce the results,
perform the following steps:

* Install R version 3.0.x or higher by following the instructions at http://www.r-project.org
* From the R environment, install and load the "devtools" package
```R
> install.packages("devtools")
> library(devtools)
```
* Install the "mixedPower" package directly from Github.com
```R
> install_github(repo="mixedPower", user="SampleSizeShop", ref="develop")
```
* Load the library
```R
> library(mixedPower)
```
* Run the simulation study (may take several hours to run with empirical calculations). You may specify output directories for data files (study.data.dir) and figures (study.figures.dir). If omitted, both default to the current working directories. To save time, you may optionally skip the empirical calculations by setting study.runEmpirical=FALSE. Note, the empirical calculations run more quickly in SAS.  If you have SAS installed on 
your local machine (and "sas.exe" on your system path), your may set study.sasEmpirical=TRUE to execute the
empirical calculations in SAS.
```R
> runSimulationStudy(study.data.dir="myDataDir", study.figures.dir="myFiguresDir", study.runEmpirical=TRUE)
```
* Run the applied example
```R
> runLongitudinalExample()
```

