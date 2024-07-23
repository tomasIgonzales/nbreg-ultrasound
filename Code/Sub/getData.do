
// Import all grayscale histograms into a single dataset

local fileList: dir "`ROOT'/Images" files "*.dta", respectcase

foreach curFile of local fileList{
	
    capture append using "`ROOT'/Images/`curFile'"
    capture use "`ROOT'/Images/`curFile'"
	
}

gen NBin_Prob = NBin/Pixels

levelsof PID,       local(PID_List)
levelsof Machine,   local(Machine_List)
levelsof Muscle,    local(Muscle_List)
