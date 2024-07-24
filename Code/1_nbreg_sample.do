version 17.0

clear
set more off
set seed 12345
log close _all
args ROOT

include "Sub/getData.do"

// Run the model in the entire dataset

log using "`ROOT'/nbreg_sample.log", name(nbreg_sample) replace text nomsg

qui{
        
    noisi di _n "Zero-inflated generalized negative binomial regression (pooled analysis)" _n

    noisi zignbreg  GSL b2.Muscle##b2.Machine [pweight=NBin_Prob]       ///
                        ,                                               ///
                        lna(b2.Muscle##b2.Machine)                      ///
                        inf(b2.Muscle##b2.Machine)                      ///
                        vce(cluster PID)             

    foreach curMuscle of local Muscle_List{

        // Define muscle labels
        if `curMuscle' == 0 local mlab = "Upper trapezius" 
        if `curMuscle' == 1 local mlab = "Pectoralis major" 
        if `curMuscle' == 2 local mlab = "Middle deltoid" 
        if `curMuscle' == 3 local mlab = "Brachioradialis" 
        if `curMuscle' == 4 local mlab = "Rectus femoris"
        if `curMuscle' == 5 local mlab = "Tibialis anterior" 

        foreach curMachine of local Machine_List{
            
            // Define machine labels
            if `curMachine' == 0 local dlab = "d0"
            if `curMachine' == 1 local dlab = "d1"
            if `curMachine' == 2 local dlab = "d2"
            
            noisi di _n "Coefficients for Machine: `dlab' Muscle: `mlab'" _n 

            // Get mu coefficients
            noisi lincom    _b[GSL:_cons]                                   +   ///
                            _b[GSL:`curMuscle'.Muscle]                      +   ///
                            _b[GSL:`curMachine'.Machine]                    +   ///
                            _b[GSL:`curMuscle'.Muscle#`curMachine'.Machine]

            local curMean_`curMachine' = exp(r(estimate))

            // Get alpha coefficients
            noisi lincom    _b[lnalpha:_cons]                                   +   /// 
                            _b[lnalpha:`curMuscle'.Muscle]                      +   ///
                            _b[lnalpha:`curMachine'.Machine]                    +   ///
                            _b[lnalpha:`curMuscle'.Muscle#`curMachine'.Machine]

            // Compute k, alpha, and p for plotting
            local curK_`curMachine' = exp(-r(estimate))
            local curA_`curMachine' = exp(r(estimate))
            local curP_`curMachine' = `curMean_`curMachine''/(`curMean_`curMachine''+`curK_`curMachine'')

            local meanLab_`curMachine'  = "`=trim("`: display %10.0f `curMean_`curMachine'''")'"
            local alphaLab_`curMachine' = "`=trim("`: display %10.2f `curA_`curMachine'''")'"

            // Get pi coefficients
            noisi lincom    _b[inflate:_cons]                                   +   /// 
                            _b[inflate:`curMuscle'.Muscle]                      +   ///
                            _b[inflate:`curMachine'.Machine]                    +   ///
                            _b[inflate:`curMuscle'.Muscle#`curMachine'.Machine]

            local curZero_`curMachine' = invlogit(r(estimate))

        }

        // Begin plotting probability density curvies
        preserve
        clear

        set obs 256
        gen GSL = _n-1
        
        foreach curMachine of local Machine_List{

            // Construct curve for current machine-muscle pair
            gen Prob_`curMachine' = .    

            forvalues i = 1/255{
            replace Prob_`curMachine' = (nbinomial(`i'      ,`curK_`curMachine'',`curP_`curMachine'') - ///
                                         nbinomial(`=`i'+1' ,`curK_`curMachine'',`curP_`curMachine'')   ///
                                        )*(1-`curZero_`curMachine'')                                    ///
                                        in `=`i'+1'
            }
        }
        
        // Since only machine 2 has zeros, we ammend label if <0.001
        local curZero_2 = "`=trim("`: display %10.3f `curZero_2''")'"
        if `curZero_2' < 0.001 local curZero_2 = "<0.001"

        // Construct plot of probability density curves for current muscle across machines
        set graphics off

        #delimit ;
        twoway  (area Prob_0 GSL, fcolor(navy%50)   lcolor(navy%0)) 
                (area Prob_1 GSL, fcolor(maroon%50) lcolor(maroon%0)) 
                (area Prob_2 GSL, fcolor(teal%50)   lcolor(teal%0)) 
                (
                ,
                ytitle("Density", size(3)) 
                xtitle("x"      , size(3)) 

                xlab(0 128 255      , nogrid labsize(3) angle(0)) 
                ylab(0 0.015 0.03   , nogrid labsize(3) angle(0)) 

                title("`mlab'", color(black) size(3))

                text(0.0290  50  "d0: {&mu}=`meanLab_0', {&alpha}=`alphaLab_0'" , size(3) color(navy)    placement(e))
                text(0.0255  50  "d1: {&mu}=`meanLab_1', {&alpha}=`alphaLab_1'" , size(3) color(maroon)  placement(e))
                text(0.0220  50  "d2: {&mu}=`meanLab_2', {&alpha}=`alphaLab_2'" , size(3) color(teal)    placement(e))
                text(0.0185  50  "      {&pi}: `curZero_2'"                     , size(3) color(teal)    placement(e))

                graphregion(color(white)) 
                legend(off)

                name(g`curMuscle', replace)
                )
                ;
        #delimit cr
        restore

    }

    log close nbreg_sample
}

// Combine into single panel and export

set graphics on

graph combine g0 g1 g2 g3 g4 g5, graphregion(color(white) margin(l=15 r=15 t=20))
graph export "`ROOT'/nbreg_sample.png" , height(2000) width(2750) replace
