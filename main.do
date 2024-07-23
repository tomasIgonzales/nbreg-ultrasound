clear
set more off

args ROOT

do nbreg_sample.do      "`ROOT'"
do nbreg_individual.do  "`ROOT'"
do nbreg_comparison.do  "`ROOT'"
do nbreg_conversion.do  "`ROOT'"
do nbreg_LOOCV.do       "`ROOT'"

