clear
set more off

args ROOT

// Install Stata user packages needed for analyses
cap net install st0336_1.pkg
cap net install qreg2.pkg

// Perform analyses
do 1_nbreg_sample.do      "`ROOT'"
do 2_nbreg_individual.do  "`ROOT'"
do 3_nbreg_comparison.do  "`ROOT'"
do 4_nbreg_conversion.do  "`ROOT'"
do 5_nbreg_LOOCV.do       "`ROOT'"

