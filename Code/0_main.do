version 17.0

/********************************************************************/
/*                                                                  */
/* GRAYSCALE HISTOGRAM HARMONIZATION ACROSS ULTRASOUND DEVICES      */
/*                                                                  */
/* This code executes a series of analyses to develop and           */
/* evaluate conversion models for grayscale histograms derived      */
/* from musculoskeletal ultrasound images across different devices. */
/*                                                                  */
/* Usage: do main_analysis.do "<ROOT_DIRECTORY>"                    */
/*                                                                  */
/* Author: Tomas Isaac Gonzales                                     */
/* Date: 24/07/2024                                                 */
/*                                                                  */
/********************************************************************/

clear
set more off
set seed 12345
args ROOT

// Install necessary Stata user packages (if not already installed)
cap net install st0336_1.pkg
cap net install qreg2.pkg



// -------------------//
// HISTOGRAM MODELING //
// -------------------//

// 1_nbreg_sample.do: 
// Zero-inflated generalized negative binomial regression to assess overall differences in 
// grayscale histogram shapes across device-muscle combinations (pooled analysis).
do 1_nbreg_sample.do        "`ROOT'"   

// 2_nbreg_individual.do:
// Zero-inflated generalized negative binomial regression to quantify variation in 
// grayscale histogram shapes (μ, α, π) across device-muscle combinations (individual-level analysis).
do 2_nbreg_individual.do    "`ROOT'"   



// --------------------------------------//
// CONVERSION MODEL DEVELOPMENT MODELING //
// --------------------------------------//

// 3_nbreg_comparison.do: 
// Sensitivity analyses to determine the optimal set of predictors for the 
// conversion models (ln μ, ln α, muscle site, interactions).
do 3_nbreg_comparison.do  "`ROOT'"

// 4_nbreg_conversion.do:
// Develop heteroskedastic linear regression models to predict ln μ and ln α for one device based on another device's values.
do 4_nbreg_conversion.do  "`ROOT'"

// 5_nbreg_LOOCV.do:
// Perform leave-one-out cross-validation (LOOCV) to assess the predictive accuracy of the conversion models and 
// visualize agreement using scatterplots, Bland-Altman plots, and box plots.
do 5_nbreg_LOOCV.do       "`ROOT'"
