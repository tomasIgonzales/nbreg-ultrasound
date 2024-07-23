clear
set more off
set seed 12345
args ROOT

include "Sub/getData.do"

qui foreach Var in lnmean lnalpha probzero{
    foreach subVar in hat lb ub{
        gen `Var'_`subVar'_gnbreg = .
    }
}

gen converged = .

// Estimate negative binomial regression model for each machine and muscle pair within each participant

qui foreach curPID of local PID_List{

    noisi di "PID: `curPID'"

    // Run model using current base level for muscle (starting at 0).  
    // Note that one participant (PID 14,) this will generate initial vals for ML estimation that cause the model to not converge.
    // In this instance, we try again using the next muscle as the base level.
    // The code below is constructed to try all base levels, but in all but one participant the first base level results in convergence.

    local baseMuscle = 0
    local hasConverged = 0

    while `hasConverged' == 0{

        #delimit ;
        noisi   zignbreg GSL b`baseMuscle'.Muscle##b2.Machine [pweight=NBin_Prob] 
                if 
                PID == `curPID'
                , 
                lna(b`baseMuscle'.Muscle##b2.Machine) 
                inf(b`baseMuscle'.Muscle##b2.Machine) 
                vce(cluster Image) 
                iterate(100)
                ;
        #delimit cr

        if e(converged) == 1    local hasConverged = 1
        else                    local baseMuscle = `baseMuscle'+1
        if `baseMuscle' == 6    local hasConverged = -1 

    }

    replace converged = e(converged) if PID == `curPID'

    // Write model coefficients

    foreach curMachine of local Machine_List{
        foreach curMuscle of local Muscle_List{
            
            // Write coefficients for mu
            lincom  _b[GSL:_cons]                                   +   ///
                    _b[GSL:`curMuscle'.Muscle]                      +   ///
                    _b[GSL:`curMachine'.Machine]                    +   ///
                    _b[GSL:`curMuscle'.Muscle#`curMachine'.Machine]
            
            replace lnmean_hat_gnbreg  =  r(estimate)   if PID == `curPID' & Machine == `curMachine' & Muscle == `curMuscle'
            replace lnmean_lb_gnbreg   =  r(lb)         if PID == `curPID' & Machine == `curMachine' & Muscle == `curMuscle'
            replace lnmean_ub_gnbreg   =  r(ub)         if PID == `curPID' & Machine == `curMachine' & Muscle == `curMuscle'

            // Write coefficients for alpha
            lincom  _b[lnalpha:_cons]                                   +   ///
                    _b[lnalpha:`curMuscle'.Muscle]                      +   ///
                    _b[lnalpha:`curMachine'.Machine]                    +   ///
                    _b[lnalpha:`curMuscle'.Muscle#`curMachine'.Machine]

            replace lnalpha_hat_gnbreg  =  r(estimate)   if PID == `curPID' & Machine == `curMachine' & Muscle == `curMuscle'
            replace lnalpha_lb_gnbreg   =  r(lb)         if PID == `curPID' & Machine == `curMachine' & Muscle == `curMuscle'
            replace lnalpha_ub_gnbreg   =  r(ub)         if PID == `curPID' & Machine == `curMachine' & Muscle == `curMuscle'

            // Write coefficients for pi
            lincom  _b[inflate:_cons]                                   +   /// 
                    _b[inflate:`curMuscle'.Muscle]                      +   ///
                    _b[inflate:`curMachine'.Machine]                    +   ///
                    _b[inflate:`curMuscle'.Muscle#`curMachine'.Machine]

            // Note that this sometimes will fail due to zero counts being very low in zome participants. We therefore try to write anyways!
            capture replace probzero_hat_gnbreg = invlogit(r(estimate)) if PID == `curPID' & Machine == `curMachine' & Muscle == `curMuscle'
            capture replace probzero_lb_gnbreg  = invlogit(r(lb))       if PID == `curPID' & Machine == `curMachine' & Muscle == `curMuscle'
            capture replace probzero_ub_gnbreg  = invlogit(r(ub))       if PID == `curPID' & Machine == `curMachine' & Muscle == `curMuscle'
          
        }  
    }
}

// Note that the model converged for all participants

keep    PID Machine Muscle lnmean* lnalpha* probzero*
order   PID Machine Muscle lnmean* lnalpha* probzero*
sort    PID Machine Muscle
duplicates drop
compress

save "`ROOT'/nbreg_individual_estimates.dta", replace
