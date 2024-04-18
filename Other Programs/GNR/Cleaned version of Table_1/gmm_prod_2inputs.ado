program gmm_prod_2inputs
	version 11

	syntax varlist [if], at(name) rhs(varlist)
	
	local m1: word 1 of `varlist'
	local m2: word 2 of `varlist'

	local vg: word 1 of `rhs'
	local k: word 2 of `rhs'
	local k2: word 3 of `rhs'
	local vg_1: word 4 of `rhs'
	local k_1: word 5 of `rhs'
	local k2_1: word 6 of `rhs'

	tempname ak ak2
	scalar `ak'=`at'[1,1]
	scalar `ak2'=`at'[1,2]

	tempvar w w_1 w2_1 w3_1 csi

	quietly gen double `w'=`vg'-`ak'*`k'-`ak2'*`k2' `if'
	quietly gen double `w_1'=`vg_1'-`ak'*`k_1'-`ak2'*`k2_1' `if'
	quietly gen double `w2_1'=`w_1'*`w_1' `if'
	quietly gen double `w3_1'=`w2_1'*`w_1' `if'
	quietly reg `w' `w_1' `w2_1' `w3_1' `if'

	quietly predict `csi' `if', resid
	
	quietly replace `m1'=`k'*`csi' `if'	
	quietly replace `m2'=`k2'*`csi' `if'	

end


