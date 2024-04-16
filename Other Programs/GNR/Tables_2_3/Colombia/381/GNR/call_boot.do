/*   GENERATE BOOTSTRAP RESULTS   */

set mem 800m
set more off
local boot_rep 200

capture log close
set logtype text
log using boot.log, replace
log off

local i=0
while `i' < `boot_rep'  {
   local i = `i' + 1
   log on
   di `i'
   log off
   do testing_boot
      if test ~= 1  {
      local i = `i' - 1
	  log on
	  di "stuck"
	  log off
	  }
   else  {
      drop test
      save boot_results_indiv_`i', replace
      clear
	  }
   }
  
   
   log close
   
   
 /*  COMPILE BOOTSTRAP RESULTS  */
   
 forvalues i=1(1)`boot_rep'  {
    use boot_results_indiv_`i'
	gen iter = `i'
	if `i' == 1  {
	   save boot_results, replace
	   }
	else  {
	   append using boot_results
	   save boot_results, replace
	   }
	erase boot_results_indiv_`i'.dta
	}
	
/*  COMPUTE STANDARD ERRORS  */

use boot_results
drop iter
foreach var of varlist *  {
   egen sd_`var' =sd(`var')
   drop `var'
   rename sd_`var' `var'
   }
collapse *, fast
save boot_results_se, replace

do testing_boot_point_est
clear
use point_est
append using boot_results_se
append using boot_results_se
foreach var of varlist   *   {
   replace `var' = `var'[1] / `var'[2] if _n == 3
   }
save boot_results_all, replace




