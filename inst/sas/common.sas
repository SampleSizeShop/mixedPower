/*
* Paper: Kenward and Roger Power
* Author: Sarah Kreidler
* Date: 4/19/2013
*
* Common.sas
*
* Set up directories for input and output.  Assumes that
* the SAS working directory is the "prog" directory.
*
*/

* directory containing modules and macros;
%LET MODULES_DIR = .;

* output directory for data sets, which are produced during data analysis;
libname outData '..\..\data';
%LET OUT_DATA_DIR = ..\..\data;

