version 17.0

clear
set more off
set seed 12345
log close _all
args ROOT

include "Sub/prepareEstimates.do"

log using "`ROOT'/nbreg_comparison.log", name(nbreg_comparison) replace text nomsg

qui{

    noisi di _n "Sensitivity analyses to determine the optimal set of predictors for the conversion models (ln μ, ln a, muscle site, interactions)." _n

    foreach Out of local Machine_List{

        // Iterate through each machine-machine combo

        local exclude = `Out'
        local others : list Machine_List - exclude

        local s1: word 1 of `others'
        local s2: word 2 of `others'

        foreach Pred in `s1' `s2'{
            
            // Iterate through outcome vars lnmean and lnalpha

            foreach Var in lnmean lnalpha{

                noisi di _n "Predict `Var' `Out' using `Pred'" _n

                // Define model levels depending on whether mu or alpha is the dependent var

                if "`Var'" == "lnmean"{
                    
                    local modelLevel_1 c.lnmean_hat_gnbreg_m`Pred'
                    local modelLevel_2 c.dummyMuscles_2 c.dummyMuscles_3 c.dummyMuscles_4 c.dummyMuscles_5 c.dummyMuscles_6 
                    local modelLevel_3 c.lnalpha_hat_gnbreg_m`Pred'
                    local modelLevel_4 c.lnmean_hat_gnbreg_m`Pred'#c.lnalpha_hat_gnbreg_m`Pred'
                }

                if "`Var'" == "lnalpha"{
                    
                    local modelLevel_1 c.lnalpha_hat_gnbreg_m`Pred'
                    local modelLevel_2 c.dummyMuscles_2 c.dummyMuscles_3 c.dummyMuscles_4 c.dummyMuscles_5 c.dummyMuscles_6
                    local modelLevel_3 c.lnmean_hat_gnbreg_m`Pred'
                    local modelLevel_4 c.lnmean_hat_gnbreg_m`Pred'#c.lnalpha_hat_gnbreg_m`Pred'
                }

                // Progressively build the model for each level

                local curModel

                forvalues i = 1/4{

                    local curModel `curModel' `modelLevel_`i''
                    noisi hetregress `Var'_hat_gnbreg_m`Out' `curModel', het(`curModel') vce(cluster PID)
                    
                    // Get log-likelihood, AIC, and whether added parameter at the current level was statistically significant

                    local logLike`i' = e(ll)
                    local AIC`i' = 2*(e(k))-2*e(ll)

                    test [`Var'_hat_gnbreg_m`Out']: `modelLevel_`i''
                    local pvalParam`i' =  r(p)

                    noisi di _n "Level `i'" _n "logLike: `logLike`i''" _n "AIC: `AIC`i''" _n "pvalParam: `pvalParam`i''" _n

                    // Output results to excel file

                    putexcel set "`ROOT'/nbreg_comparison", sheet("ln`Var'_m`Out' using m`Pred'") modify

                    local j = (`i'-1)*4

                    putexcel A`=1+`j'' = ("Model_`i'")

                    putexcel A`=2+`j'' = ("logLike")
                    putexcel A`=3+`j'' = ("AIC")
                    putexcel A`=4+`j'' = ("pvalParam")

                    putexcel B`=2+`j'' = ("`=trim("`: display %10.2f `logLike`i'''")'")
                    putexcel B`=3+`j'' = ("`=trim("`: display %10.2f `AIC`i'''")'")
                    putexcel B`=4+`j'' = ("`=cond(`pvalParam`i''<0.001,"<0.001","`=trim("`: display %10.3f `pvalParam`i'''")'")'")
            
                }     
            }
        }
    }

    log close nbreg_comparison

}
