version 17.0

// Import grayscale histogram data file into Stata

import delimited using "`ROOT'/Data.csv", case(preserve) asdouble

// Reshape into long format

keep PID Machine Muscle Image Pixels NBin*
egen tempID = concat(PID Machine Muscle Image), punct("_")
reshape long NBin, i(tempID) j(GSL)
drop tempID

sort    PID Machine Muscle Image GSL
order   PID Machine Muscle Image Pixels GSL NBin

gen NBin_Prob = NBin/Pixels

// Initialize locals for different levels of ID, Machine, and Muscle

levelsof PID,       local(PID_List)
levelsof Machine,   local(Machine_List)
levelsof Muscle,    local(Muscle_List)
