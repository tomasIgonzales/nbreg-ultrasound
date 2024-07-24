version 17.0

use "`ROOT'/nbreg_individual_estimates.dta", clear

// Reshape the individual esitmates dataset into wide format for conversion modelling

keep PID Machine Muscle *_hat_*

rename *_hat_* =_m
egen id2 = concat(PID Muscle), punct(_)
reshape wide *_hat_*, i(id2) j(Machine)

drop id2 
order PID Muscle
sort  PID Muscle

local Machine_List 0 1 2
levelsof Muscle, local(Muscle_List)
levelsof PID, local(PID_List)

tab Muscle, gen(dummyMuscles_)

foreach Var in lnmean lnalpha{
    foreach curMachine of local Machine_List{
        egen `Var'_std_gnbreg_m`curMachine' = std(`Var'_hat_gnbreg_m`curMachine')
    }
}


