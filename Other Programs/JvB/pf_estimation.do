cd "${gwp}/programs/3_production_function_estimation"

clear
set more off

*1-data for variable values
local nbs		=	"${gwp}/intermediate_data/NBSworking1998-2007.dta"
local merged	=	"${gwp}/intermediate_data/merged.dta"
*take STA values from merged.dta, i.e the maxim in STA has been capped by NBS obs.

*2-sample choice (density ratio specification and random seed in dr estiamtion for simulated samples)
local sample	=	"${gsample}"

if "`sample'"~="STAoriginal" & "`sample'"~="NBSoriginal" {
	local drseed	=	$gdrseed
	local kn		=	$gkn
	local sampledir	=	"${gwp}/intermediate_data/simulated_samples/`sample'_seed_`drseed'_kn`kn'"
}

*3-years of samples for estimation: 1st stage and 2nd stage moment condition specification
local mc		=	"1yr"
*local mc		=	"2yrs"

if "`mc'"=="1yr" & "`sample'"=="NBSoriginal" {	
	local sample1st	=	""
	local lag		=	1
	local sample2nd	=	""
	local gmmprog	=	"GNR_gmm"
}

if "`mc'"=="1yr" & "`sample'"~="NBSoriginal" {	
	if "${gkeep2008}" == "no" {
		local sample1st	=	"keep if time>=2009"
		local lag		=	1
		local sample2nd	=	"for var vg l k ll kk lk: replace  X_1 = . if time==2009"
		*2009-2013 sample will be used in the 1st stage fitting of materials share
		*since 2008 is problematic, lag values of 2009 not reliable, thus not used in the 2nd stage; replace with missing
		*2010,2011,2012 and 2013 will be used in the 2nd stage; note 2007 is not used either
	}
	if "${gkeep2008}" == "yes" {
		local sample1st	=	""
		local lag		=	1
		local sample2nd	=	""	
	}	
	
	local gmmprog	=	"GNR_gmm"
	*note that the evolution of omega includes only up to the 3rd order of the lag
}

/*
if "`mc'"=="2yrs" & "`sample'"~="NBSoriginal" {	
	local sample1st	=	"keep if time~=2008"
	local lag		=	2
	local sample2nd	=	""
	*2007 and 2009-2013 sample will be used in the 1st stage fitting of materials share
	*with two periods lag and exclusion of 2008, the lag values for 2010 are automatically missing; do nothing
	*2009, 2011, 2012 and 2013 will be used in the second stage

	local gmmprog	=	"GNR_gmm_22"
	*note that the evolution of omega includes only up to the 2nd order of the lag
	*with higher orders it takes longer and is less likely to achieve convergence
}
*/

*4-naming files for saving coef estimation results
if "`sample'"=="NBSoriginal" | "`sample'"=="STAoriginal" {
	local sim0		=	0
	local sim1		=	0
	local coefresult=	"${gwp}/intermediate_results/pf_estimation/coef_`sample'_mc_`mc'_lag.dta"
	local estresult=	"${gwp}/intermediate_results/pf_estimation/estimates/estimates_`sample'_mc_`mc'_lag.dta"
}
else {
	local sim0		=	$gs0
	local sim1		=	$gs1
	local coefresult=	"${gwp}/intermediate_results/pf_estimation/coef_sim`sim0'-`sim1'_`sample'_seed_`drseed'_kn_`kn'_mc_`mc'_lag.dta"
	local estresult=	"${gwp}/intermediate_results/pf_estimation/estimates/estimates_sim`sim0'-`sim1'_`sample'_seed_`drseed'_kn_`kn'_mc_`mc'_lag.dta"
}

*5-input deflators
local deflator "MovingMS"

*****************************************************************************************

*Prepare variable values for production function estimation

if "`sample'"=="NBSoriginal" use "`nbs'",clear
else  {
	use TAX* year if TAXoutput>=5000 & TAXoutput~=. & int(TAXcic_adj/100)>=13 & int(TAXcic_adj/100)<=42 & int(TAXcic_adj/100)~=16 & int(TAXcic_adj/100)~=28 using "`merged'", clear 
	rename TAX* *
	rename id TAXid
	gen input=output-va+va_tax
	gen imt`deflator'=log(input/inputdeflator`deflator')
	gen l=log(employment)
	gen k=log(real_cap)

	gen ownership="1SOE" if type=="110" | type=="141" | type=="143" | type=="151" | substr(type,1,2)=="16"
	replace ownership="3HMT" if substr(type,1,1)=="2"
	replace ownership="4FOR" if substr(type,1,1)=="3"
	replace ownership="2NonSOE" if ownership=="" 
	
	gen province=substr(dq,1,2)
	gen cic2=int(cic_adj/100)
	
	duplicates drop TAXid year, force
}

quietly{
	for var imt`deflator' \ any i : rename X Y
	gen mshare=(input/output)
	replace mshare=2 if mshare>2 & mshare~=.
	gen si = ln(mshare)
	gen yg = log(output/outputdeflator)
	
	cap gen id=TAXid
	cap gen id=NBSid
	
	keep id year yg si mshare i l k ownership province cic2 cic_adj
	keep if (si~=. & si<3) & l~=. & k~=. & i~=. & yg~=. 
	
	gen ll =l*l
	gen kk =k*k
	gen ii =i*i
	gen lk =l*k
	gen li =l*i
	gen ki =k*i
	gen lki=l*k*i
} /*quiely implemented*/
rename id PFid
local pfvariables:  tempfile
save `pfvariables', replace

*Estimate production function

forval s=`sim0'(1)`sim1' {

	if "`sample'"~="STAoriginal" & "`sample'"~="NBSoriginal" {
	
		di "iteration " `s'

		use year TAXid using "`sampledir'/sample_`sample'_`s'.dta",clear	
		*use the original TAXid rather than the id in the dr estimation input files
		
		quietly {
			bysort TAXid year: gen duplicate=_n
			
			gen sampled=1
			reshape wide sampled, i(TAXid duplicate) j(year)
			reshape long sampled, i(TAXid duplicate) j(year)
			bysort TAXid duplicate (year): keep if sampled==1|sampled[_n+1]==1
			
			gen PFid=TAXid
			
			merge n:1 year PFid using `pfvariables', keep(matched) nogenerate
			bysort PFid duplicate (year): keep if sampled==1|sampled[_n+1]==1  /* drops 2013 obs. if only 2014 was sampled */
			
			egen   id = group(PFid duplicate)
			rename year time 
			drop duplicate 
		}/*quiely implemented*/
		
		local temp_one_sim: tempfile
		save `temp_one_sim'
	}/*for simulated sample; skipped with original*/

	forval cc=13(1)41 {
		di "  cic2=  `cc'"
		if `cc'==16 | `cc'==28 | `cc'==38 continue
		
		if "`sample'"=="STAoriginal" | "`sample'"=="NBSoriginal"  {
			use if cic2==`cc' using `pfvariables', clear
			rename year time
			egen id=group(PFid)
		}

		else use if cic2==`cc' using `temp_one_sim', clear

		quietly {
			
			di "*** 1st step ***"
			`sample1st'
			
			for var si mshare \ any -3 0.05: replace X=Y if X<Y
			regress mshare l k i ll kk ii lk li ki lki
			predict crap
			replace crap = crap - _b[_cons]
			egen mcrap = min(crap)
			scalar ncrap=mcrap
			drop crap mcrap
			scalar ncrap=-ncrap + 0.01
			
			capture nl ( si = ln( {g0=ncrap} + {gl=_b[l]}*l + {gk=_b[k]}*k + {gi=_b[i]}*i + {gll=_b[ll]}*ll + {glk=_b[lk]}*lk + ///
					 {gli=_b[li]}*li + {gkk=_b[kk]}*kk + {gki=_b[ki]}*ki + {gii=_b[ii]}*ii + {glki=_b[lki]}*lki) ), iter(200) delta(1e-10)
			if e(converge)~= 1 exit
			
			predict melas 
			predict eg, resid
			replace eg=-eg
			egen mexp_eg=mean(exp(eg))
			replace melas=melas-ln(mexp_eg)
			replace melas=exp(melas)
			mat beta=e(b)
			svmat double beta
			for num 1/11 \ any 0 l k i ll lk li kk ki ii lki: rename betaX gY
			for            any 0 l k i ll lk li kk ki ii lki: egen    GX = mean(gX)
			for            any 0 l k i ll lk li kk ki ii lki: replace GX = GX / mexp_eg
			clear matrix
			
			gen integ_G_I = i * ( G0+Gl*l+Gk*k + Gll*ll+Gkk*kk+Glk*lk + (Gi*i+Gli*li+Gki*ki+Glki*lki)/2 + Gii*ii/3)
			gen vg        = yg - eg - integ_G_I
			
			*** 2nd step ***
			
			tabulate ownership, gen(OO)
			tabulate time     , gen(YY)
			tabulate province , gen(PP)
			tabulate cic_adj  , gen(CC)
			
			for var vg l k ll kk lk: regress X OO* YY* PP* CC* \ predict PX, res \ summ X \ replace PX = PX + _result(3)
			tset id time
			for var vg l k ll kk lk: gen  X_`lag' = L`lag'.PX
			`sample2nd'
			
			reg Pvg Pl Pk Pll Pkk Plk 
			matrix test2 = e(b)
			matrix test3 = test2[1,1..5]
			matrix test3[1,1] = test2[1,"Pl"]
			matrix test3[1,2] = test2[1,"Pk"]
			matrix test3[1,3] = test2[1,"Pll"]
			matrix test3[1,4] = test2[1,"Pkk"]
			matrix test3[1,5] = test2[1,"Plk"]
			matrix drop test2
			
			gmm `gmmprog' if Pvg!=. & vg_`lag'!=. & l_`lag'!=. & k_`lag'!=., one nequations(5) parameters(al ak all akk alk) from(test3) winitial(identity) rhs(Pvg Pl Pk Pll Pkk Plk vg_`lag' l_`lag' k_`lag' ll_`lag' kk_`lag' lk_`lag') conv_maxiter(100)
			
			mat beta=e(b)
			svmat double beta
			for num 1/5 \ any l k ll kk lk: rename betaX aY
			for           any l k ll kk lk: egen   AX=mean(aX)
			
			quietly for var yg vg i l k ll kk lk: regress X OO* /*YY*/ PP* CC* \ predict QX, res
			gen logomega=Qvg-Al*Ql-Ak*Qk-All*Qll-Akk*Qkk-Alk*Qlk
			gen tfp     = vg-Al* l-Ak* k-All* ll-Akk* kk-Alk* lk
			tset id time
			gen tfpg  = logomega-L.logomega
			gen tfpg2 = tfp -L.tfp
			
			gen lelas = Gl*i+2*Gll*li+Glk*ki + Gli*ii/2 + Glki*ki*i/2 + Al + 2*All*l + Alk*k
			gen kelas = Gk*i+2*Gkk*ki+Glk*li + Gki*ii/2 + Glki*li*i/2 + Ak + 2*Akk*k + Alk*l
			gen rts   = melas + lelas + kelas
			gen sim=`s'
						
			preserve 
				for var i yg \ any m q: rename X Y
				keep id time cic_adj cic2 ownership province mshare q l k m lelas kelas melas rts tfp tfpg* sim
			
				if `cc'~=13 | `s'~=`sim0' append using  "`estresult'"
				save "`estresult'", replace
			restore
						
			keep sim cic2 G* A* melas lelas kelas rts tfp tfpg*
			collapse (median) G* A* m=melas l=lelas k=kelas rts tfp tfpg* (mean) Mm=melas Ml=lelas Mk=kelas Mrts=rts Mtfp=tfp Mtfpg=tfpg Mtfpg2=tfpg2 ///
					 (sd) sdm=melas sdl=lelas sdk=kelas sdrts=rts sdtfp=tfp sdtfpg=tfpg sdtfpg2=tfpg2, by(cic2 sim)					 
			
		}/*quietly implemented lines*/

		if `cc'~=13 | `s'~=`sim0' append using "`coefresult'"
		save "`coefresult'",replace
	}/*cic2 loop*/
}/*sim loop*/
