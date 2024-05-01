/************************************************************/
/* First step: Estimate share regression with NLLS 			*/
/************************************************************/
/**  The code below runs the share regression using a      **/ 
/**  a (log)polynomial approximation.  This step     	   **/
/**  recovers the output elasticity of flexible inputs 	   **/
/**                                                        **/
/**  The initial values are set using an OLS regression    **/
/**  of shares (si) on a polynomial in inputs.  If the     **/
/**  non-linear least squares procedure (nl) fails to      **/
/**  converge, the initial values can be changed.	   	   **/
/**********************************************************/ {
use "C:\Users\marku\Documents\GNRProdEst\Other Programs\GNR\Cleaned version of Table_1\cd_data_500.dta", clear
* Get starting values from linear regression
cap drop crap
regress si k i kk ii ki if si~=. & k~=. & i~=. // Linear regression of what we want to get non-linearly
predict crap if si~=. & k~=. & i~=.  // Get fitted values
replace crap = crap - _b[_cons] // Substract constant
cap drop mcrap 
egen mcrap = min(crap) // Take the minimum
scalar ncrap=mcrap
drop crap mcrap
scalar ncrap=-ncrap + 0.1 // Take negative of the minimum and add 0.1. No idea why but ok

scalar list ncrap

#delimit;
nl ( si = ln({g0=ncrap} + {gk=_b[k]}*k + {gi=_b[i]}*i + {gkk=_b[kk]}*kk + {gki=_b[ki]}*ki + {gii=_b[ii]}*ii)) if si~=. & k~=. & i~=. ;
#delimit cr

cap drop ielas eg
predict ielas if k~=. & si~=. & i~=.
predict eg if k~=. & si~=. & i~=., resid

replace eg=-eg
replace ielas=exp(ielas)
* The next part only puts the coefficients into the dataset as variables
mat beta=e(b)
svmat double beta
ren beta1 g0
ren beta2 gk
ren beta3 gi
ren beta4 gkk
ren beta5 gki
ren beta6 gii

foreach var of varlist g0-gii {
	egen s`var'=mean(`var')
	drop `var'
	ren s`var' `var'
}	
clear matrix

* Calculate the integral now
gen integ_G_I = g0+gk*k+gkk*kk + (gi*i+gki*ki)/2 + gii*ii/3 // They never correct the coefficients. This is different than in the paper and in R-version of the command. I do not think this is correct.
replace integ_G_I=integ_G_I*i // Because in equation it says m^(r+1) but we did not add the + 1 in the part above
gen vg = yg - eg - integ_G_I
exit
tset id time
gen vg_1=L.vg
gen k_1=L.k
gen kk_1=L.kk
}
/************************************************************/
/**  Now to recover the remaining coefficients associated  **/
/**  with capital and labor (the constant of the PDE)      **/
/**                                                        **/
/**  The initial values are set using an OLS regression    **/
/**  of vg on log capital and loglabor                     **/
/************************************************************/

reg vg k kk 
matrix test2 = e(b)
matrix test3 = test2[1,1..2]
matrix test3[1,1] = test2[1,"k"]
matrix test3[1,2] = test2[1,"kk"]
matrix drop test2

gmm gmm_prod_2inputs if vg!=. & k!=. & vg_1!=. & k_1!=., one nequations(2) parameters (ak akk) from(test3) winitial(identity) rhs(vg k kk vg_1 k_1 kk_1)
mat beta=e(b)
svmat double beta
ren beta1 ak
ren beta2 akk
foreach var of varlist ak-akk {
	egen s`var'=mean(`var')
	drop `var'
	ren s`var' `var'
}	
clear matrix

/************************************************************/
/**  Generate productivity in logs and levels, and compute **/
/**  the output elasticities of labor and capital.         **/
/************************************************************/

gen logomega=vg-ak*k-akk*kk
gen omega=exp(logomega)

gen kelas=gk*i+2*gkk*ki + gki*ii/2 + ak + 2*akk*k

/************************************************************/
/**  End of the code.                                      **/
/************************************************************/

