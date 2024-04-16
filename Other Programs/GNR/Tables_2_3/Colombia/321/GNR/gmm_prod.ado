program gmm_prod
	version 11

	syntax varlist [if], at(name) rhs(varlist)
	
	local m1: word 1 of `varlist'
	local m2: word 2 of `varlist'
	local m3: word 3 of `varlist'
	local m4: word 4 of `varlist'
	local m5: word 5 of `varlist'

	local vg: word 1 of `rhs'
	local l: word 2 of `rhs'
	local k: word 3 of `rhs'
	local l2: word 4 of `rhs'
	local k2: word 5 of `rhs'
	local lk: word 6 of `rhs'
	local vg_1: word 7 of `rhs'
	local l_1: word 8 of `rhs'
	local k_1: word 9 of `rhs'
	local l2_1: word 10 of `rhs'
	local k2_1: word 11 of `rhs'
	local lk_1: word 12 of `rhs'

	tempname al ak al2 ak2 alk
	scalar `al'=`at'[1,1] 
	scalar `ak'=`at'[1,2]
	scalar `al2'=`at'[1,3]
	scalar `ak2'=`at'[1,4]
	scalar `alk'=`at'[1,5]

	tempvar w w_1 w2_1 w3_1 csi

	quietly gen double `w'=`vg'-`al'*`l'-`ak'*`k'-`al2'*`l2'-`ak2'*`k2'-`alk'*`lk' `if'
	quietly gen double `w_1'=`vg_1'-`al'*`l_1'-`ak'*`k_1'-`al2'*`l2_1'-`ak2'*`k2_1'-`alk'*`lk_1' `if'
	quietly gen double `w2_1'=`w_1'*`w_1' `if'
	quietly gen double `w3_1'=`w2_1'*`w_1' `if'
	quietly reg `w' `w_1' `w2_1' `w3_1' `if'
	quietly predict `csi' `if', resid
	
	quietly replace `m1'=`l'*`csi' `if'	
	quietly replace `m2'=`k'*`csi' `if'	
	quietly replace `m3'=`l2'*`csi' `if'	
	quietly replace `m4'=`k2'*`csi' `if'	
	quietly replace `m5'=`lk'*`csi' `if'	

end

