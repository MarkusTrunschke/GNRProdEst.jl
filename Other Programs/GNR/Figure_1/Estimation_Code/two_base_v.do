clear
 
local periods "3 5 10 20 30 50"
local firms "200 500"
foreach per of local periods {
foreach fir of local firms  {
use ..\Monte_Carlo_Data\two_base_v_data, clear
rename yt prod

keep if t <=`per'+150
keep if i <= `fir'



matrix ini1 = (0.65,0.25)

gen test_m_1 = .



quietly  {

forvalues i=1(1)500  {

preserve
keep if s == `i'

rename i id
rename t time
tsset id time

gen y = ln(prod)
gen k = ln(kt)
gen l = ln(mt)

gen l2 = l*l
gen k2 = k*k
gen lk = l*k


reg y lrp k l k2 l2 lk
predict eps, resid

gen double vg = y - eps

gen double vg_1 = l.vg
gen double l_1 = l.l
gen double k_1 = l.k


gen al_m_1 = .
gen ak_m_1 = .
capture gmm gmm_prod_proxy if vg!=. & l!=. & k!=. & vg_1!=. & l_1!=. & k_1!=., one nequations(2) parameters(al ak) from(ini1) winitial(identity) rhs(vg l k vg_1 l_1 k_1) conv_maxiter(50)
replace test_m_1 = e(converged)
replace al_m_1 = _b[/al]
replace ak_m_1 = _b[/ak]



keep if _n == 1
keep s al* ak* test*

if `i'==1  {
save two_base_v_mc_`fir'_`per', replace
}
else  {
append using two_base_v_mc_`fir'_`per'
save two_base_v_mc_`fir'_`per', replace
}

restore
}

}
}




}
