*LOOP OVER CD, TL (translog), and CES SIMULATIONS


*CES VERSION
clear all
set more off


clear all
set obs 1
gen ss = 0

save boot2_ces, replace

local sims 100
forvalues j=1(1)`sims'  {

infile id time h kt mt yt wt using ..\Monte_Carlo_Data\ces_data.out, clear
keep if id <= 500*`j' & id > 500*(`j'-1)


replace time=time-170
gen si_level=mt/yt
ren yt yg_level
ren kt k_level
ren mt i_level
gen lkt=log(k_level)
gen lmt=log(i_level)
gen lyt=log(yg_level)
gen lwt=log(wt)
gen lk2=lkt*lkt
gen lm2=lmt*lmt
gen lkm=lkt*lmt

gen arg=0.25*((exp(lkt))^0.5)+0.65*((exp(lmt))^0.5)
gen kelas_true=0.9*0.25*((exp(lkt))^0.5) / arg
gen ielas_true=0.9*0.65*((exp(lmt))^0.5) / arg

do GNR_code_2inputs

egen mkelas=mean(kelas)
egen mielas=mean(ielas)
egen sdkelas=sd(kelas)
egen sdielas=sd(ielas)

egen mkelas_t=mean(kelas_t)
egen mielas_t=mean(ielas_t)
egen sdkelas_t=sd(kelas_t)
egen sdielas_t=sd(ielas_t)

gen sumelas = kelas + ielas
gen sumelas_t = kelas_t + ielas_t
egen msumelas=mean(sumelas)
egen msumelas_t=mean(sumelas_t)
egen sdsumelas=sd(sumelas)
egen sdsumelas_t=sd(sumelas_t)


gen tempbkelas = kelas< 0 | kelas > 1
gen tempbielas = ielas< 0 | ielas > 1

gen tempbkelas_t = kelas_t< 0 | kelas_t > 1
gen tempbielas_t = ielas_t< 0 | ielas_t > 1

egen bkelas = mean(tempbkelas) 
egen bielas = mean(tempbielas) 
egen bkelas_t = mean(tempbkelas_t) 
egen bielas_t = mean(tempbielas_t) 


keep mkelas* mielas* sdkelas* sdielas* bkelas* bielas* msumelas* sdsumelas* 
keep if _n == 1

	  gen ss = `j'
	  sleep 100
	  append using boot2_ces
	  sleep 100
      save boot2_ces, replace
      clear

}
//////////////////////////////////////



*TL VERSION
clear all
set more off


clear all
set obs 1
gen ss = 0

save boot2_tl, replace


local sims 100
forvalues j=1(1)`sims'  {

infile id time h kt mt yt wt using ..\Monte_Carlo_Data\tl_data.out, clear
keep if id <= 500*`j' & id > 500*(`j'-1)


replace time=time-170
gen si_level=mt/yt
ren yt yg_level
ren kt k_level
ren mt i_level
gen lkt=log(k_level)
gen lmt=log(i_level)
gen lyt=log(yg_level)
gen lwt=log(wt)
gen lk2=lkt*lkt
gen lm2=lmt*lmt
gen lkm=lkt*lmt

gen kelas_true=0.25+0.015*2*lkt - 0.032*lmt
gen ielas_true=0.65+0.015*2*lmt - 0.032*lkt


do GNR_code_2inputs


egen mkelas=mean(kelas)
egen mielas=mean(ielas)
egen sdkelas=sd(kelas)
egen sdielas=sd(ielas)

egen mkelas_t=mean(kelas_t)
egen mielas_t=mean(ielas_t)
egen sdkelas_t=sd(kelas_t)
egen sdielas_t=sd(ielas_t)

gen sumelas = kelas + ielas
gen sumelas_t = kelas_t + ielas_t
egen msumelas=mean(sumelas)
egen msumelas_t=mean(sumelas_t)
egen sdsumelas=sd(sumelas)
egen sdsumelas_t=sd(sumelas_t)


gen tempbkelas = kelas< 0 | kelas > 1
gen tempbielas = ielas< 0 | ielas > 1

gen tempbkelas_t = kelas_t< 0 | kelas_t > 1
gen tempbielas_t = ielas_t< 0 | ielas_t > 1

egen bkelas = mean(tempbkelas) 
egen bielas = mean(tempbielas) 
egen bkelas_t = mean(tempbkelas_t) 
egen bielas_t = mean(tempbielas_t) 


keep mkelas* mielas* sdkelas* sdielas* bkelas* bielas* msumelas* sdsumelas* 
keep if _n == 1

	  gen ss = `j'
	  sleep 100
	  append using boot2_tl
	  sleep 100
      save boot2_tl, replace
      clear

}
/////////////////////////////////////


*CD VERSION
clear all
set more off


clear all
set obs 1
gen ss = 0

save boot2_cd, replace


local sims 100
forvalues j=1(1)`sims'  {



infile id time h kt mt yt using ..\Monte_Carlo_Data\cd_data.out, clear
keep if id <= 500*`j' & id > 500*(`j'-1)

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


do GNR_code_2inputs

egen mkelas=mean(kelas)
egen mielas=mean(ielas)
egen sdkelas=sd(kelas)
egen sdielas=sd(ielas)

egen mkelas_t=mean(kelas_t)
egen mielas_t=mean(ielas_t)
egen sdkelas_t=sd(kelas_t)
egen sdielas_t=sd(ielas_t)

gen sumelas = kelas + ielas
gen sumelas_t = kelas_t + ielas_t
egen msumelas=mean(sumelas)
egen msumelas_t=mean(sumelas_t)
egen sdsumelas=sd(sumelas)
egen sdsumelas_t=sd(sumelas_t)


gen tempbkelas = kelas< 0 | kelas > 1
gen tempbielas = ielas< 0 | ielas > 1

gen tempbkelas_t = kelas_t< 0 | kelas_t > 1
gen tempbielas_t = ielas_t< 0 | ielas_t > 1

egen bkelas = mean(tempbkelas) 
egen bielas = mean(tempbielas) 
egen bkelas_t = mean(tempbkelas_t) 
egen bielas_t = mean(tempbielas_t) 


keep mkelas* mielas* sdkelas* sdielas* bkelas* bielas* msumelas* sdsumelas*
keep if _n == 1

	  gen ss = `j'
	  sleep 100
	  append using boot2_cd
	  sleep 100
      save boot2_cd, replace
      clear

}
/////////////////////////////////////





