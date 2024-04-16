clear
local periods "3 5 10 20 30 50"
local firms "200 500"
local varis "half_base_ base_ two_base_ ten_base_"
foreach vari of local varis  { 
	do `vari'v
	foreach per of local periods  {
		foreach fir of local firms  {
			use `vari'v_mc_`fir'_`per', clear  
			quietly replace al_m_1 = . if test_m_1 == 0
			quietly replace ak_m_1 = . if test_m_1 == 0		
			foreach var of varlist al_m_1 ak_m_1   {
				quietly centile `var', centile(2.5 97.5)
				quietly gen double ci_l_`var'= r(c_1)
				quietly gen double ci_h_`var'= r(c_2)
				quietly egen m_`var' = mean(`var')
			}  
			quietly egen m_test_m_1 = sum(test_m_1)
			gen firms = `fir'
			gen periods = `per'
			gen str20 variation = "`vari'"
			keep firms periods variation ci* m_*
			keep if _n == 1
			if `per' == 3 & `fir'== 200 & "`vari'" == "half_base_"  {
				save allv_results, replace
			}
			else  {
				append using allv_results
				sleep 3000
				save allv_results, replace
			}
		}
	}
}		




clear
use allv_results
gen var_num = .
replace var_num = 1 if variation == "half_base_"
replace var_num = 2 if variation == "base_"
replace var_num = 3 if variation == "two_base_"
replace var_num = 4 if variation == "ten_base_"
label define levels 1 "Half Base" 2 "Base" 3 "Two Base" 4 "Ten Base" 		
label values var_num levels
		
sort var_num firms periods			

order variation firms periods m_test_m_1 m_al_m_1 ci_l_al_m_1 ci_h_al_m_1
keep variation firms periods m_test_m_1 m_al_m_1 ci_l_al_m_1 ci_h_al_m_1
li variation firms periods m_test_m_1 m_al_m_1 ci_l_al_m_1 ci_h_al_m_1, ab(20)

outsheet using mc_time_series, comma replace
