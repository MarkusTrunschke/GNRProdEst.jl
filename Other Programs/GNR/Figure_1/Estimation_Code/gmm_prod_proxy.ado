program gmm_prod_proxy
	version 11

	syntax varlist [if], at(name) rhs(varlist)
	
	local m1: word 1 of `varlist'
	local m2: word 2 of `varlist'

	local vg: word 1 of `rhs'
	local l: word 2 of `rhs'
	local k: word 3 of `rhs'
	local vg_1: word 4 of `rhs'
	local l_1: word 5 of `rhs'
	local k_1: word 6 of `rhs'

	tempname al ak
	scalar `al'=`at'[1,1] 
	scalar `ak'=`at'[1,2]

	tempvar w w_1 csi 

	quietly gen double `w'=`vg'-`al'*`l'-`ak'*`k' `if'
	quietly gen double `w_1'=`vg_1'-`al'*`l_1'-`ak'*`k_1' `if'
	
	quietly reg `w' `w_1' `if'
	quietly predict `csi' `if', resid 
	
	quietly replace `m1'=`l_1'*`csi' `if'	
	quietly replace `m2'=`k'*`csi' `if'	

	
end

