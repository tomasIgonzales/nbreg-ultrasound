version 17.0

clear
set more off
set seed 12345
log close _all
args ROOT

include "Sub/prepareEstimates.do"

log using "`ROOT'/4_nbreg_conversion.log", name(nbreg_conversion) replace text nomsg

qui{

    noisi di _n "Develop heteroskedastic linear regression models to predict ln μ and ln α for one device based on another device's values." _n

    foreach Out of local Machine_List{

        // Iterate through each machine-machine combo

        local exclude = `Out'
        local others : list Machine_List - exclude

        local s1: word 1 of `others'
        local s2: word 2 of `others'

        foreach Pred in `s1' `s2'{
            
            qui foreach Var in lnmean lnalpha{

                noisi di _n "Predict `Var' `Out' using `Pred'" _n
                
                /////////////////////////////////////////////////
                // Run model for current machine-machine combo //
                /////////////////////////////////////////////////

                noisi hetregress    `Var'_hat_gnbreg_m`Out'             ///
                                    c.lnmean_hat_gnbreg_m`Pred'         ///  
                                    c.lnalpha_hat_gnbreg_m`Pred'        ///
                                    b4.Muscle                           ///
                                    ,                                   ///
                                    het(                                ///
                                        c.lnmean_std_gnbreg_m`Pred'     ///
                                        c.lnalpha_std_gnbreg_m`Pred'    /// 
                                        b4.Muscle                       ///
                                        )                               ///
                                    vce(cluster PID)

                // Compute spearman rank correlation and RMSE for model
                    
                predict `Var'_pred

                qui spearman `Var'_pred `Var'_hat_gnbreg_m`Out' 
                local `Var'_rho = "`=trim("`: display %10.3f r(rho)'")'" 

                gen  `Var'_sqerror = (`Var'_pred - `Var'_hat_gnbreg_m`Out')^2
                qui su `Var'_sqerror
                local `Var'_rmse = "`=trim("`: display %10.3f sqrt(r(mean))'")'"
                
                drop `Var'_pred `Var'_sqerror

                noisi di _n "Spearman: ``Var'_rho'" _n "RMSE: ``Var'_rmse'" _n

                /////////////////////////////////////////////////
                // Prepare excel spreadsheet for model outputs //
                /////////////////////////////////////////////////

                putexcel set "`ROOT'/nbreg_conversion", sheet("`Var'_m`Out' using m`Pred'") modify
            
                putexcel A1 = ("Outcome")
                putexcel A2 = ("`Var' from m`Out'")

                putexcel B1 = ("Using")
                putexcel B2 = ("Predictors from m`Pred'")

                // Write spearman rank correlation
                putexcel C1 = ("Spearman")
                putexcel C2 = ("`=trim("`: display %10.3f ``Var'_rho''")'")

                // Write RMSE
                putexcel D1 = ("RMSE_model")
                putexcel D2 = ("`=trim("`: display %10.3f ``Var'_rmse''")'")

                // Write column headers
                putexcel A4 = ("predictor")
                putexcel B4 = ("estimate")
                putexcel C4 = ("se")
                putexcel D4 = ("lb_95CI")
                putexcel E4 = ("ub_95CI")
                putexcel F4 = ("p_value")
                putexcel G4 = ("Formatted")


                //////////////////////////////////////////
                // Write coefficients for mean function //
                //////////////////////////////////////////

                noisi di _n  "Write coefficients for mean function" _n
                
                // Write coefficient for mu
                noisi lincom _b[`Var'_hat_gnbreg_m`Out': c.lnmean_hat_gnbreg_m`Pred']

                putexcel A5 = ("lnmean_m`Pred'")
                putexcel B5 = ("`=trim("`: display %10.3f r(estimate)'")'")
                putexcel C5 = ("`=trim("`: display %10.3f r(se)'")'")
                putexcel D5 = ("`=trim("`: display %10.3f r(lb)'")'")
                putexcel E5 = ("`=trim("`: display %10.3f r(ub)'")'")
                putexcel F5 = ("`=trim("`: display %10.3f r(p)'")'")
                putexcel G5 = (`"`=trim("`: display %10.2f r(estimate)'")' (`=trim("`: display %10.2f r(se)'")')"')
                
                // Write coefficient for alpha
                noisi lincom _b[`Var'_hat_gnbreg_m`Out': c.lnalpha_hat_gnbreg_m`Pred']

                putexcel A6 = ("lnalpha_m`Pred'")
                putexcel B6 = ("`=trim("`: display %10.3f r(estimate)'")'")
                putexcel C6 = ("`=trim("`: display %10.3f r(se)'")'")
                putexcel D6 = ("`=trim("`: display %10.3f r(lb)'")'")
                putexcel E6 = ("`=trim("`: display %10.3f r(ub)'")'")
                putexcel F6 = ("`=trim("`: display %10.3f r(p)'")'")
                putexcel G6 = (`"`=trim("`: display %10.2f r(estimate)'")' (`=trim("`: display %10.2f r(se)'")')"')


                // Write coefficient for constant intercept
                #delimit ;
                noisi lincom  (
                                (_b[`Var'_hat_gnbreg_m`Out': _cons] + _b[`Var'_hat_gnbreg_m`Out': 0.Muscle]) +
                                (_b[`Var'_hat_gnbreg_m`Out': _cons] + _b[`Var'_hat_gnbreg_m`Out': 1.Muscle]) +
                                (_b[`Var'_hat_gnbreg_m`Out': _cons] + _b[`Var'_hat_gnbreg_m`Out': 2.Muscle]) +
                                (_b[`Var'_hat_gnbreg_m`Out': _cons] + _b[`Var'_hat_gnbreg_m`Out': 3.Muscle]) +
                                (_b[`Var'_hat_gnbreg_m`Out': _cons] + _b[`Var'_hat_gnbreg_m`Out': 4.Muscle]) +
                                (_b[`Var'_hat_gnbreg_m`Out': _cons] + _b[`Var'_hat_gnbreg_m`Out': 5.Muscle]) 
                              ) / 6
                ;
                #delimit cr
                
                putexcel A7 = ("Intercept")
                putexcel B7 = ("`=trim("`: display %10.3f r(estimate)'")'")
                putexcel C7 = ("`=trim("`: display %10.3f r(se)'")'")
                putexcel D7 = ("`=trim("`: display %10.3f r(lb)'")'")
                putexcel E7 = ("`=trim("`: display %10.3f r(ub)'")'")
                putexcel F7 = ("`=trim("`: display %10.3f r(p)'")'")
                putexcel G7 = (`"`=trim("`: display %10.2f r(estimate)'")' (`=trim("`: display %10.2f r(se)'")')"')

                // Write coefficients for muscle-specific intercepts
                foreach curMuscle of local Muscle_List{

                    #delimit ;
                    noisi lincom    _b[`Var'_hat_gnbreg_m`Out': _cons]              + 
                                    _b[`Var'_hat_gnbreg_m`Out': `curMuscle'.Muscle] - 
                                    (
                                        (
                                            (_b[`Var'_hat_gnbreg_m`Out': _cons] + _b[`Var'_hat_gnbreg_m`Out': 0.Muscle]) +
                                            (_b[`Var'_hat_gnbreg_m`Out': _cons] + _b[`Var'_hat_gnbreg_m`Out': 1.Muscle]) +
                                            (_b[`Var'_hat_gnbreg_m`Out': _cons] + _b[`Var'_hat_gnbreg_m`Out': 2.Muscle]) +
                                            (_b[`Var'_hat_gnbreg_m`Out': _cons] + _b[`Var'_hat_gnbreg_m`Out': 3.Muscle]) +
                                            (_b[`Var'_hat_gnbreg_m`Out': _cons] + _b[`Var'_hat_gnbreg_m`Out': 4.Muscle]) +
                                            (_b[`Var'_hat_gnbreg_m`Out': _cons] + _b[`Var'_hat_gnbreg_m`Out': 5.Muscle]) 
                                        )
                                        / 6 
                                    )
                    ;
                    #delimit cr

                    putexcel A`=9+`curMuscle'' = ("`curMuscle'.Muscle")
                    putexcel B`=9+`curMuscle'' = ("`=trim("`: display %10.3f r(estimate)'")'")
                    putexcel C`=9+`curMuscle'' = ("`=trim("`: display %10.3f r(se)'")'")
                    putexcel D`=9+`curMuscle'' = ("`=trim("`: display %10.3f r(lb)'")'")
                    putexcel E`=9+`curMuscle'' = ("`=trim("`: display %10.3f r(ub)'")'")
                    putexcel F`=9+`curMuscle'' = ("`=trim("`: display %10.3f r(p)'")'")
                    putexcel G`=9+`curMuscle'' = (`"`=trim("`: display %10.2f r(estimate)'")' (`=trim("`: display %10.2f r(se)'")')"')

                }


                //////////////////////////////////////////////
                // Write coefficients for variance function //
                //////////////////////////////////////////////

                noisi di _n  "Write coefficients for variance function" _n

                putexcel A16 = ("Variance function chi2 p_value")
                putexcel A17 = ("`=trim("`: display %10.3f e(p_c)'")'")

                // Write coefficient for mu
                noisi lincom (_b[lnsigma2: c.lnmean_std_gnbreg_m`Pred'])

                putexcel A20 = ("lnmean_m`Pred'")
                putexcel B20 = ("`=trim("`: display %10.2f r(estimate)'")'")
                putexcel C20 = ("`=trim("`: display %10.2f r(se)'")'")
                putexcel D20 = ("`=trim("`: display %10.2f r(lb)'")'")
                putexcel E20 = ("`=trim("`: display %10.2f r(ub)'")'")
                putexcel F20 = ("`=trim("`: display %10.2f r(p)'")'")

                // Write coefficient for alpha
                noisi lincom (_b[lnsigma2: c.lnalpha_std_gnbreg_m`Pred']) 

                putexcel A21 = ("lnalpha_m`Pred'")
                putexcel B21 = ("`=trim("`: display %10.2f r(estimate)'")'")
                putexcel C21 = ("`=trim("`: display %10.2f r(se)'")'")
                putexcel D21 = ("`=trim("`: display %10.2f r(lb)'")'")
                putexcel E21 = ("`=trim("`: display %10.2f r(ub)'")'")
                putexcel F21 = ("`=trim("`: display %10.2f r(p)'")'")

                // Write coefficient for constant intercept
                #delimit ;
                noisi lincom  (
                                (_b[lnsigma2: _cons] + _b[lnsigma2: 0.Muscle]) +
                                (_b[lnsigma2: _cons] + _b[lnsigma2: 1.Muscle]) +
                                (_b[lnsigma2: _cons] + _b[lnsigma2: 2.Muscle]) +
                                (_b[lnsigma2: _cons] + _b[lnsigma2: 3.Muscle]) +
                                (_b[lnsigma2: _cons] + _b[lnsigma2: 4.Muscle]) +
                                (_b[lnsigma2: _cons] + _b[lnsigma2: 5.Muscle]) 
                              ) / 6
                ;
                #delimit cr
                
                putexcel A22 = ("Intercept")
                putexcel B22 = ("`=trim("`: display %10.3f r(estimate)'")'")
                putexcel C22 = ("`=trim("`: display %10.3f r(se)'")'")
                putexcel D22 = ("`=trim("`: display %10.3f r(lb)'")'")
                putexcel E22 = ("`=trim("`: display %10.3f r(ub)'")'")
                putexcel F22 = ("`=trim("`: display %10.3f r(p)'")'")

                
                // Write coefficients for muscle-specific intercepts
                foreach curMuscle of local Muscle_List{

                    #delimit ;
                    noisi lincom    _b[lnsigma2: _cons]              + 
                                    _b[lnsigma2: `curMuscle'.Muscle] - 
                                    (
                                        (
                                            (_b[lnsigma2: _cons] + _b[lnsigma2: 0.Muscle]) +
                                            (_b[lnsigma2: _cons] + _b[lnsigma2: 1.Muscle]) +
                                            (_b[lnsigma2: _cons] + _b[lnsigma2: 2.Muscle]) +
                                            (_b[lnsigma2: _cons] + _b[lnsigma2: 3.Muscle]) +
                                            (_b[lnsigma2: _cons] + _b[lnsigma2: 4.Muscle]) +
                                            (_b[lnsigma2: _cons] + _b[lnsigma2: 5.Muscle]) 
                                        )
                                        / 6 
                                    )
                    ;
                    #delimit cr

                    putexcel A`=24+`curMuscle'' = ("`curMuscle'.Muscle")
                    putexcel B`=24+`curMuscle'' = ("`=trim("`: display %10.2f r(estimate)'")'")
                    putexcel C`=24+`curMuscle'' = ("`=trim("`: display %10.2f r(se)'")'")
                    putexcel D`=24+`curMuscle'' = ("`=trim("`: display %10.2f r(lb)'")'")
                    putexcel E`=24+`curMuscle'' = ("`=trim("`: display %10.2f r(ub)'")'")
                    putexcel F`=24+`curMuscle'' = ("`=trim("`: display %10.2f r(p)'")'")

                }
            }
        }
    }

    log close nbreg_conversion

}

