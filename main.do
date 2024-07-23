clear
set more off

args ROOT

// nbreg_sample.do
// Apply zero-inflated generali  
do nbreg_sample.do      "`ROOT'"


do nbreg_individual.do  "`ROOT'"
do nbreg_comparison.do  "`ROOT'"
do nbreg_conversion.do  "`ROOT'"
do nbreg_LOOCV.do       "`ROOT'"

