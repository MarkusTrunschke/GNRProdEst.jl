clear all
set mem 1000m
set more off

use ../data_col

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
 
keep logRGO logRVA si logRI logK logL export import id year adv hiwag 
sort id year

ren logRVA yv
ren logRGO yg
ren logRI i
ren logK k
ren logL l
foreach v of varlist * {
	replace `v'=. if `v'<-999
}
gen ll=l*l
gen kk=k*k
gen ii=i*i
gen lk=l*k
gen li=l*i
gen ki=k*i
gen lll=ll*l
gen kkk=kk*k
gen iii=ii*i
gen lki=l*k*i



reg yg l k i ll kk ii lk li ki

gen rflelasgo=_b[l] + 2*_b[ll]*l + _b[lk]*k + _b[li]*i 
egen c_mlelas_go=mean(rflelasgo)
gen rfkelasgo=_b[k] + 2*_b[kk]*k + _b[lk]*l + _b[ki]*i
egen c_mkelas_go=mean(rfkelasgo)
gen rfielasgo=_b[i] + 2*_b[ii]*i + _b[li]*l + _b[ki]*k
egen c_mielas_go=mean(rfielasgo)
gen rfsumgo=rflelasgo+rfkelasgo+rfielasgo
egen c_msum_go=mean(rfsumgo)
gen c_mKelasLelas_go=c_mkelas_go/c_mlelas_go
drop rflelasgo rfkelasgo rfielasgo rfsumgo 

predict rflprodgo, resid
gen rfprodgo=exp(rflprodgo)


centile rfprodgo, cen(10)
scalar den2=r(c_1)
gen crapgo=rfprodgo/den2

qreg crapgo, q(90)

    gen c_ratio_90_10_go = _b[_cons]


drop crapgo


centile rfprodgo, cen(5)
scalar den2=r(c_1)
gen crapgo=rfprodgo/den2

qreg crapgo, q(95)

    gen c_ratio_95_5_go = _b[_cons]

drop crapgo



centile rfprodgo, cen(25)
scalar den2=r(c_1)
gen crapgo=rfprodgo/den2

qreg crapgo, q(75)

    gen c_ratio_75_25_go = _b[_cons]

drop crapgo


sum rfprodgo if export==0
scalar den2=r(mean)
gen crapgo=rfprodgo/den2

reg crapgo export

    gen c_export_go = _b[export]

drop crapgo 


sum rfprodgo if import==0
scalar den2=r(mean)
gen crapgo=rfprodgo/den2

reg crapgo import

    gen c_import_go = _b[import]


drop crapgo




sum rfprodgo if adv==0
scalar den2=r(mean)
gen crapgo=rfprodgo/den2

reg crapgo adv

    gen c_adv_go = _b[adv]

drop crapgo



sum rfprodgo if hiwag==0
scalar den2=r(mean)
gen crapgo=rfprodgo/den2

reg crapgo hiwag

    gen c_hiwag_go = _b[hiwag]

drop crapgo


    keep c_* 
    keep if _n == 1
