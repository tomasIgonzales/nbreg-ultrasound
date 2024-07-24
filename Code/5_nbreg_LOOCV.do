version 17.0

clear
set more off
set seed 12345
log close _all
args ROOT

include "Sub/prepareEstimates.do"

log using "`ROOT'/5_nbreg_LOOCV.log", name(nbreg_LOOCV) replace text nomsg

qui{

    noisi di _n "Perform leave-one-out cross-validation (LOOCV) to assess the predictive accuracy of the conversion models" _n

    foreach Out of local Machine_List{

        // Iterate through each machine-machine combo

        local exclude = `Out'
        local others : list Machine_List - exclude

        local s1: word 1 of `others'
        local s2: word 2 of `others'

        foreach Pred in `s1' `s2'{
            
            foreach curPID of local PID_List{
            
                foreach Var in lnmean lnalpha{
                    
                    // Run model excluding one participant, then predict in that participant

                    qui hetregress     `Var'_hat_gnbreg_m`Out'              ///
                                        c.lnmean_hat_gnbreg_m`Pred'         /// 
                                        c.lnalpha_hat_gnbreg_m`Pred'        ///
                                        b4.Muscle                           ///
                                        if                                  ///
                                        PID != `curPID'                     ///
                                        ,                                   ///
                                        het(                                ///
                                            c.lnmean_std_gnbreg_m`Pred'     /// 
                                            c.lnalpha_std_gnbreg_m`Pred'    /// 
                                            b4.Muscle                       ///
                                            )                               ///
                                        vce(cluster PID)

                    predict `Var'_pred_`curPID' if PID == `curPID', equation(`Var'_hat_gnbreg_m`Out')
            
                }
            }


            foreach Var in lnmean lnalpha{
                
                noisi di _n "Predict `Var' `Out' using `Pred'" _n

                // Compute RMSE of pred vs measured

                egen `Var'_pred = rowtotal(`Var'_pred_*)  // Collapse predictions from all folds of LOOCV
                gen  `Var'_sqerror = (`Var'_pred - `Var'_hat_gnbreg_m`Out')^2

                su `Var'_sqerror
                local RMSE_loocv = "`=trim("`: display %10.3f sqrt(r(mean))'")'"

                // Convert units to their real values, and then compute spearman rank corr

                gen outReal = exp(`Var'_hat_gnbreg_m`Out')
                gen predReal = exp(`Var'_pred)

                spearman outReal predReal
                local corrcoef = "`=trim("`: display %10.2f r(rho)'")'" 

                noisi di _n "Spearman: `corrcoef'" _n "RMSE: `RMSE_loocv'" _n


                // Define ranges for labels on plots, depending on outcome variables

                if "`Var'" == "lnmean"{

                    if `Out' == 0 local labs 20 90 160
                    if `Out' == 1 local labs 50 150 250
                    if `Out' == 2 local labs 0 40 80

                    if `Out' == 0 local range 20 160
                    if `Out' == 1 local range 50 250
                    if `Out' == 2 local range 0 80

                    if `Out' == 0 local lim 160
                    if `Out' == 1 local lim 250
                    if `Out' == 2 local lim 80

                    local baLabs -1.0 -0.5 0 0.5 1.0

                    local labOut ="`=ustrunescape("\u0075")'"
                    local labPred ="`=ustrunescape("\u00FB")'"

                }

                if "`Var'" == "lnalpha"{

                    if `Out' == 0 local labs 0 0.2 0.4
                    if `Out' == 1 local labs 0 0.1 0.2
                    if `Out' == 2 local labs 0 1.5 3

                    if `Out' == 0 local range 0 0.4
                    if `Out' == 1 local range 0 0.2
                    if `Out' == 2 local range 0 3

                    if `Out' == 0 local lim 0.4
                    if `Out' == 1 local lim 0.2
                    if `Out' == 2 local lim 3

                    local baLabs -1.5 -0.75 0 0.75 1.5

                    //Note that lim above only removes plotting of one outlier for device 2

                    local labOut ="`=ustrunescape("\u0061")'"
                    local labPred ="`=ustrunescape("\u00E2")'"

                }

                // Prepare scatter plot demonstrating agreement

                set graphics off
                #delimit ;

                twoway
                (
                scatter outReal
                        predReal
                        if
                        predReal < `lim' & outReal < `lim'
                        ,
                        mlcolor(navy%0) 
                        mcolor(navy%40) 
                        msize(0.7) 
                )
                (
                func y=x
                , 
                lcolor(black%40) 
                range(`range') 
                )
                (
                ,

                title("`labOut' for d`Out' using d`Pred'", size(3) color(black))

                xtitle("`labPred'", size(3) orient(horiz))
                ytitle("`labOut'", size(3) orient(horiz))

                xlabels(`labs', nogrid labsize(3) angle(0)) 
                ylabels(`labs', nogrid labsize(3) angle(0))
                note("RMSLE: `RMSE_loocv'" "{&rho}: `corrcoef'", ring(0) position(11))
                graphregion(color(white)) 
                legend(off) 
                name(sc_`Var'_`Out'`Pred', replace)
                )
                ;
                #delimit cr


                // Prepare Bland-Altman plot

                gen diffReal = ln(predReal/outReal)
                gen aveReal = (predReal + outReal)/2

                su diffReal

                local curMean = ("`=trim("`: display %10.3f r(mean)'")'")
                local curSD = ("`=trim("`: display %10.3f r(sd)'")'")
                local uSD = r(mean) + 1.96*r(sd)
                local lSD = r(mean) - 1.96*r(sd)

                #delimit ;

                scatter diffReal 
                        aveReal
                        if
                        predReal < `lim' & outReal < `lim'
                        , 
                        mlcolor(navy%0) 
                        mcolor(navy%40) 
                        msize(0.7) 

                        yline(`curMean', lcolor(black%40)) 
                        yline(`uSD', lcolor(black%40))  
                        yline(`lSD', lcolor(black%40))

                        title("`labOut' for d`Out' using d`Pred'", size(3) color(black))

                        xtitle("Mean of `labPred' and `labOut'", size(3))
                        ytitle("ln(`labPred'/`labOut')", size(3))

                        xlabels(`labs', nogrid labsize(3) angle(0))  
                        ylabels(`baLabs', nogrid labsize(3) angle(0))
                        note("Mean: `curMean'" "SD: `curSD'", ring(0) position(1)) 
                        graphregion(color(white))
                        name(ba_`Var'_`Out'`Pred', replace)
                        ;

                #delimit cr


                // Prepare Box plots of agreement across muscle sites  
                // Check if any residuals by muscle site are statistically significantly different from zero

                noisi di _n "Assess differential agreement by muscle site" _n

                noisi qreg2 diffReal i.Muscle, cluster(PID) quantile(0.50)

                noisi lincom _b[_cons] + _b[0.Muscle]
                noisi lincom _b[_cons] + _b[1.Muscle]
                noisi lincom _b[_cons] + _b[2.Muscle]
                noisi lincom _b[_cons] + _b[3.Muscle]
                noisi lincom _b[_cons] + _b[4.Muscle]
                noisi lincom _b[_cons] + _b[5.Muscle]

                noisi lincom  ( (_b[_cons] + _b[0.Muscle]) +    ///
                                (_b[_cons] + _b[1.Muscle]) +    ///
                                (_b[_cons] + _b[2.Muscle]) +    ///
                                (_b[_cons] + _b[3.Muscle]) +    ///
                                (_b[_cons] + _b[4.Muscle]) +    ///
                                (_b[_cons] + _b[5.Muscle])      ///
                              ) / 6

                #delimit ;

                graph   box diffReal
                        if
                        predReal < `lim' & outReal < `lim'
                        ,
                        over(Muscle, label(labsize(2.5) angle(90)) relabel(1 "UT" 2 "PM" 3 "MD" 4 "BR" 5 "RF" 6 "TA"))

                        noout
                        note("")

                        title("`labOut' for d`Out' using d`Pred'", size(3) color(black))

                        yline(0,lcolor(black%40))
                        ytitle("ln(`labPred'/`labOut')", size(3))

                        ylabels(`baLabs', nogrid labsize(3) angle(0))
                        graphregion(color(white))
                        name(bx_`Var'_`Out'`Pred', replace)
                        ;

                #delimit cr
                    
                drop `Var'_pred* `Var'_sqerror outReal predReal diffReal aveReal

            }
        }
    }

    log close nbreg_LOOCV

}

set graphics on

// Combine graphs into panel and export

foreach Var in lnmean lnalpha{

    // Combine scatter plots

    #delimit ;
    graph combine   sc_`Var'_01 
                    sc_`Var'_02 
                    sc_`Var'_10 
                    sc_`Var'_12 
                    sc_`Var'_20 
                    sc_`Var'_21
                    , 
                    row(3) 
                    col(2) 
                    title("A)", position(11) just(left) size(2.5) color(black)) 
                    graphregion(color(white)) 
                    name("sc_`Var'_all", replace)
                    ;
    #delimit cr

    // Combine BA plots

    #delimit ;
    graph combine   ba_`Var'_01 
                    ba_`Var'_02 
                    ba_`Var'_10 
                    ba_`Var'_12 
                    ba_`Var'_20 
                    ba_`Var'_21
                    , 
                    row(3) 
                    col(2)
                    title("B)", position(11) just(left) size(2.5) color(black))
                    graphregion(color(white))
                    name("ba_`Var'_all", replace)
                    ;
    #delimit cr

    // Combine box plots

    #delimit ;
    graph combine   bx_`Var'_01 
                    bx_`Var'_02 
                    bx_`Var'_10 
                    bx_`Var'_12 
                    bx_`Var'_20 
                    bx_`Var'_21
                    , 
                    row(3) 
                    col(2)
                    title("C)", position(11) just(left) size(2.5) color(black))
                    graphregion(color(white))
                    name("bx_`Var'_all", replace)
                    ;
    #delimit cr


    // Combine all plots into a single panel figure
    
    #delimit ;
    graph combine   sc_`Var'_all 
                    ba_`Var'_all 
                    bx_`Var'_all
                    , 
                    col(3) 
                    imargin(0 0 0 0) 
                    graphregion(color(white) 
                    margin(l=0 r=0 t=12.5 b=12.5)) 
                    name("`Var'_all", replace) 
                    ;
    #delimit cr 

    graph export "`ROOT'/`Var'_LOOCV.png" , height(2000) width(2750) replace 

}
