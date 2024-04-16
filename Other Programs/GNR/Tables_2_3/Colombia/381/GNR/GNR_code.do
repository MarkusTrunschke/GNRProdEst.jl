

/************************************************************/
/**  This code is designed to estimate a gross output      **/
/**  production function with three inputs: capital,       **/
/**  labor, and intermediate inputs, using a series        **/
/**  estimator for the elasticity (from the share):        **/
/**  s = ln(g0+gl*l+gk*k+gi*i+gll*l*l+glk*l*k+gli*li+	   **/
/**    	 gkk*k*k+gki*k*i+gii*i*i+glki*l*k*i) + epsilon	   **/
/**  and a series estimator for the remaining part of the  **/
/**  the production function:                              **/
/**                                                        **/
/**  y= Integral(G(k,l,i)/I)dI + al*l+ak*k+all*ll+akk*kk   **/
/**    	+alk*lk+omega+epsilon                              **/
/************************************************************/


/************************************************************/
/**  Rename variables using the instructions below.        **/ 
/**  Pay attention to upper and lower case, as it          **/
/**  matters for Stata.  Inputs, output, and the share of  **/
/**  intermediate expenditures in total revenue should     **/
/**  all be expressed in levels.  The code will create     **/
/**  the log values.  Also note that input and output      **/
/**  variables are real values, whereas the share is       **/
/**  the nominal share.                                    **/
/**                                                        **/
/**  Firm ID: id                                           **/
/**  Time series variable (e.g., year, month): time        **/
/**  Real gross output: yg_level                           **/
/**  Real labor: l_level                                   **/
/**  Real capital: k_level                                 **/
/**  Real intermediate inputs): i_level                    **/
/**  Nominal Share of intermediates: si_level              **/ 
/**                                                        **/             
/************************************************************/


gen yg = ln(yg_level)
gen l = ln(l_level)
gen k = ln(k_level)
gen i = ln(i_level)
gen si = ln(si_level)

gen ll=l*l
gen kk=k*k
gen ii=i*i
gen lk=l*k
gen li=l*i
gen ki=k*i
gen lki=l*k*i



/************************************************************/
/**  The code below runs the share regression using a      **/ 
/**  a (log)polynomial approximation.  This step     	   **/
/**  recovers the output elasticity of flexible inputs 	   **/
/**                                                        **/
/**  The initial values are set using an OLS regression    **/
/**  of shares (si) on a polynomial in inputs.  If the     **/
/**  non-linear least squares procedure (nl) fails to      **/
/**  converge, the initial values can be changed.	   **/
/************************************************************/

regress si l k i ll kk ii lk li ki lki if si~=. & l~=. & k~=. & i~=.
matrix test = e(b)
predict crap if si~=. & l~=. & k~=. & i~=.
replace crap = crap - _b[_cons]
egen mcrap = min(crap)
scalar ncrap=mcrap
drop crap mcrap
scalar ncrap=-ncrap + 0.1

#delimit;
capture nl ( si = ln({g0=ncrap} + {gl=_b[l]}*l + {gk=_b[k]}*k + {gi=_b[i]}*i + {gll=_b[ll]}*ll + 
   {glk=_b[lk]}*lk + {gli=_b[li]}*li + {gkk=_b[kk]}*kk + {gki=_b[ki]}*ki +
   {gii=_b[ii]}*ii + {glki=_b[lki]}*lki) ) if si~=. & l~=. & k~=. & i~=., iter(100);
#delimit cr

gen test = e(converge)
if test ~= 1 exit

predict ielas if l~=. & k~=. & si~=. & i~=.
predict eg if l~=. & k~=. & si~=. & i~=., resid
replace eg=-eg
egen mexp_eg=mean(exp(eg))
replace ielas=ielas-ln(mexp_eg)
replace ielas=exp(ielas)
mat beta=e(b)
svmat double beta
ren beta1 g0
ren beta2 gl
ren beta3 gk
ren beta4 gi
ren beta5 gll
ren beta6 glk
ren beta7 gli
ren beta8 gkk
ren beta9 gki
ren beta10 gii
ren beta11 glki
foreach var of varlist g0-glki {
	egen s`var'=mean(`var')
	drop `var'
	ren s`var' `var'
	replace `var' = `var' / mexp_eg
}	
clear matrix

gen integ_G_I = g0+gl*l+gk*k+gll*ll+gkk*kk+glk*lk + (gi*i+gli*li+gki*ki+glki*lki)/2 + gii*ii/3
replace integ_G_I=integ_G_I*i 
gen vg = yg - eg - integ_G_I

tset id time
gen vg_1=L.vg
gen l_1=L.l
gen k_1=L.k
gen ll_1=L.ll
gen kk_1=L.kk
gen lk_1=L.lk
gen vg_2=L2.vg
gen l_2=L2.l
gen k_2=L2.k
gen ll_2=L2.ll
gen kk_2=L2.kk
gen lk_2=L2.lk

/************************************************************/
/**  Now to recover the remaining coefficients associated  **/
/**  with capital and labor (the constant of the PDE)      **/
/**                                                        **/
/**  The initial values are set using an OLS regression    **/
/**  of vg on log capital and loglabor                     **/
/************************************************************/

reg vg l k ll kk lk 
matrix test2 = e(b)
matrix test3 = test2[1,1..5]
matrix test3[1,1] = test2[1,"l"]
matrix test3[1,2] = test2[1,"k"]
matrix test3[1,3] = test2[1,"ll"]
matrix test3[1,4] = test2[1,"kk"]
matrix test3[1,5] = test2[1,"lk"]
matrix drop test2

capture gmm gmm_prod if vg!=. & l!=. & k!=. & vg_1!=. & l_1!=. & k_1!=., one nequations(5) parameters(al ak all akk alk) from(test3) winitial(identity) rhs(vg l k ll kk lk vg_1 l_1 k_1 ll_1 kk_1 lk_1) conv_maxiter(100)
replace test = e(converged)
if test ~= 1 exit


mat beta=e(b)
svmat double beta
ren beta1 al
ren beta2 ak
ren beta3 all
ren beta4 akk
ren beta5 alk
foreach var of varlist al-alk {
	egen s`var'=mean(`var')
	drop `var'
	ren s`var' `var'
}	

/************************************************************/
/**  Generate productivity in logs and levels, and compute **/
/**  the output elasticities of labor and capital.         **/
/************************************************************/

gen logomega=vg-al*l-ak*k-all*ll-akk*kk-alk*lk
gen omega=exp(logomega)

gen lelas=gl*i+2*gll*li+glk*ki + gli*ii/2 + glki*ki*i/2 + al + 2*all*l + alk*k
gen kelas=gk*i+2*gkk*ki+glk*li + gki*ii/2 + glki*li*i/2 + ak + 2*akk*k + alk*l

/************************************************************/
/**  End of the code.                                      **/
/************************************************************/
