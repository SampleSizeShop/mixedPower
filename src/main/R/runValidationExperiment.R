#
# Validation experiment for the manuscript describing
# power for the Wald test with Kenward Roger denominator
# degrees of freedom in the fixed balanced and fixed
# unbalanced case.
#
# Author: Sarah Kreidler
# Created Date: 8/11/2013
#
# !This program requires that the current working directory
# !is set to the location of runValidationExperiment.R
#
# Inputs:
# ../
#
# Outputs:
#
#
library(lme4)
library(plyr)

#
# set some relative paths
#
INPUT_DIR = "../../../input/";
OUTPUT_DATA_DIR = "../../../output/datasets/";
MODULES_DIR = "../../modules/R/";

#
# Define the studyDesign class
#
source(paste(MODULES_DIR,"studyDesign.R",sep=""))
#
# load the SAS simulation module
#
source(paste(MODULES_DIR,"sasEmpiricalPower.R"),sep="")

#
# Load the study designs
#
# This defines a list of designs in the variable
# 'studyDesignList'
#
source(paste(INPUT_DIR,"studyDesignList.R"))


#
#
# calculateApproximateKRPower
# 
# Function applied to each object in the studyDesignList
# Performs the following functions
# 1. Calls SAS from the command line to obtain empirical power.
#    Empirical power is written to a csv file in the output/datasets
#    directory
# 2. 
#