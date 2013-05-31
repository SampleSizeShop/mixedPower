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
%LET MODULES_DIR = ..\..\modules\sas;

* input directory for data sets, including raw and manually edited data;
libname inData '..\..\..\input';
%LET IN_DATA_DIR = ..\..\..\input;

* output directory for data sets, which are produced during data analysis;
libname outData '..\..\..\output\datasets';
%LET OUT_DATA_DIR = ..\..\..\output\datasets;

* output directory for figures ;
libname outFig '..\..\..\output\figures';
%LET OUT_FIGURES_DIR = ..\..\..\output\figures;

* output directory for tables ;
libname outTab '..\..\..\output\tables';
%LET OUT_TABLES_DIR = ..\..\..\output\tables;

* output directory for automated reports, logs, and listings;
libname outRep '..\..\..\output\reports';
%LET OUT_REPORT_DIR = ..\..\..\output\reports;
