*CD VERSION
clear all
set more off
pause on
/* Get a subset of the Monte Carlo data */ {
	infile id time h kt mt yt using "C:\Users\marku\Documents\GNRProdEst\Other Programs\GNR\Cleaned version of Table_1\cd_data.out", clear
	keep if id <= 2000
	replace time=time-170
	gen si_level=mt/yt
	ren yt yg_level
	ren kt k_level
	ren mt i_level
	gen lkt=log(k_level)
	gen lmt=log(i_level)
	gen lyt=log(yg_level)
	gen lk2=lkt*lkt
	gen lm2=lmt*lmt
	gen lkm=lkt*lmt
	gen kelas_true=0.25
	gen ielas_true=0.65
	gen yg = ln(yg_level)
	gen k = ln(k_level)
	gen i = ln(i_level)
	gen si = ln(si_level)
	gen kk=k*k
	gen ii=i*i
	gen ki=k*i

	save "C:\Users\marku\Documents\GNRProdEst\Other Programs\GNR\Cleaned version of Table_1\cd_data_2000.dta", replace
	export delimited "C:\Users\marku\Documents\GNRProdEst\Other Programs\GNR\Cleaned version of Table_1\cd_data_2000.csv", replace
}
use "C:\Users\marku\Documents\GNRProdEst\Other Programs\GNR\Cleaned version of Table_1\cd_data_500.dta", clear

* Run estimation
do GNR_code_2inputs_clean


* Mean is the important statistic here (comparable to ACF point estimate)
sum kelas ielas


