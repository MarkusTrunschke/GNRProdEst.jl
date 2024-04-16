clear all
set more off

use ../data_chi

/*
   program tsset_boot
      sort id year
      bys id: gen first = year[_n]==year[_n-1]
      bys id year: gen second = sum(first)
      gen double id2 = 0.01*second + id
      format id2 %8.2f
      rename id idold
      rename id2 id
	 drop first second idold
   end program

   bsample, cluster(id)
   tsset_boot
*/   

keep id year RGO L K RI si export import hiwag adv 
foreach v of varlist * {
	replace `v'=. if `v'<-999
}
replace si=exp(si)
ren year time
ren RGO yg_level
ren L l_level
ren K k_level
ren RI i_level
ren si si_level

do GNR_code

gen logprod=logomega+eg

gen prod=exp(logprod)

egen mlelas=mean(lelas)
egen mkelas=mean(kelas)
egen mielas=mean(ielas)

gen sumelas=lelas+kelas+ielas
egen msumelas=mean(sumelas)
gen mKelasmLelas=mkelas/mlelas
drop lelas kelas ielas sumelas 


sum mielas mlelas mkelas msumelas mKelasmLelas


centile prod, cen(10)
scalar den2=r(c_1)
gen crap=prod/den2

capture qreg crap, q(90)
replace test = 1 if e(convcode) == 0
replace test = . if e(convcode) ~= 0 
if test ~= 1 exit



     gen c_str_ratio_90_10 =_b[_cons]


drop crap

centile prod, cen(5)
scalar den2=r(c_1)
gen crap=prod/den2

capture qreg crap, q(95)
replace test = 1 if e(convcode) == 0
replace test = . if e(convcode) ~= 0 
if test ~= 1 exit



     gen c_str_ratio_95_5 =_b[_cons]


drop crap

centile prod, cen(25)
scalar den2=r(c_1)
gen crap=prod/den2

capture qreg crap, q(75)
replace test = 1 if e(convcode) == 0
replace test = . if e(convcode) ~= 0 
if test ~= 1 exit



     gen c_str_ratio_75_25 =_b[_cons]
	 

drop crap

sum prod if export==0
scalar den2=r(mean)
gen crap=prod/den2

reg crap export

     gen c_str_export =_b[export]

drop crap 

sum prod if import==0
scalar den2=r(mean)
gen crap=prod/den2

reg crap import

     gen c_str_import =_b[import]

drop crap 



sum prod if adv==0
scalar den2=r(mean)
gen crap=prod/den2

reg crap adv

     gen c_str_adv =_b[adv]

drop crap 

sum prod if hiwag==0
scalar den2=r(mean)
gen crap=prod/den2

reg crap hiwag

     gen c_str_hiwag =_b[hiwag]

drop crap 



     keep test mlelas mkelas mielas msumelas mKelas c_*
     keep if _n == 1
	 
	 save point_est, replace
	 
	 
	 
	 
	 
	 
	 
