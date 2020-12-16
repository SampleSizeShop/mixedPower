The mixedPower Package
=========================

The mixedPower package for R (>4.0.0) calculates power for the 
linear mixed model for longitudinal and cluster randomized designs. 
The power calculations may be used for designs with or without 
anticipated missing data.  

The package provides companion code for the manuscript:

Kreidler, S. M., Ringham, B. M., Muller, K. E., & Glueck, D. H. 
A Power Approximation for the Kenward and Roger Wald Test 
in the Linear Mixed Mode, In review.

### Instructions for replicating the manuscript results 

The results in the above manuscript were produced using R version 4.0.2. 
To reproduce the results, perform the following steps:

* Install R version 4.0.x or higher by following the instructions at http://www.r-project.org
* From the R environment, install and load the "devtools" package
```R
> install.packages("devtools")
> library(devtools)
```
* Install the "rPowerlib" package from Github.com
```R
> install_github("SampleSizeShop/rPowerlib")
```
* Install the "invWishartSum" package from Github.com
```R
> install_github("SampleSizeShop/invWishartSum")
```
* Install the "mixedPower" package from Github.com
```R
> install_github("SampleSizeShop/mixedPower")
```
* Load the library
```R
> library(mixedPower)
```
* Run the simulation study. By default, data files and figures will be written to the
current working directory. The study.data.dir and study.figures.dir arguments can 
be used to override the defaults. Empirical results require a SAS (v9.4) installation.
They take several hours to run, and so are turned off by default. A copy of the empirical 
results are included in the R package.
```R
> runSimulationStudy(study.data.dir="myDataDir", study.figures.dir="myFiguresDir", study.runEmpirical=FALSE)
```
* Run the applied example
```R
> runLongitudinalExample(study.data.dir="myDataDir", study.figures.dir="myFiguresDir")
```

