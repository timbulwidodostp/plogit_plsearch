* Penalised logistic regression (lasso)
program define plogi_l2
	version 8.2
	args todo b lnf g negH
	tempvar theta lnfj d2
	tempname bsq penal sb llnp P1 P2 delta s el bj
	mleval `theta'=`b'
	quietly gen double `lnfj'=-ln(1+exp(-`theta')) if $ML_y1==1
	quietly replace `lnfj'=-ln(1+exp(`theta')) if $ML_y1==0
	mlsum `lnf'=`lnfj'
	scalar `llnp' = `lnf'
	* Penalisation
	local c=colsof(`b')
	local j 1
	* Delta controls region in which GA's funky fudge comes into play
	scalar `penal'=0
	while `j'<`c' {
        scalar `bj'=`b'[1,`j']
		scalar `s'=sd[1,`j']
*		if abs(`sb'[1,`j'])>$delta {
			scalar `penal'=`penal'+$lambda*abs(`bj'*`s')
*		}
*		else {
*			scalar `penal'=`penal'+$lambda*(3/8*$delta+6*(`sb'[1,`j'])^2/(8*$delta)-(`sb'[1,`j'])^4/(8*$delta^3))
*		}
		local j=`j'+1
	}
	scalar `lnf'=`lnf'-`penal'
	* !! PR bug fix
/*
	return scalar llnp = `llnp'
	return scalar llpen = `penal'
*/
	global PL_llnp=`llnp'
	global PL_llpen=`penal'
	if `todo'==0 | `lnf'==. exit

	tempvar d1
	* First derivatives
	quietly gen double `d1'=1/(exp(`theta')+1) if $ML_y1==1
	quietly replace `d1'=1/(exp(`theta')+1)-1 if $ML_y1==0
	mlvecsum `lnf' `g' = `d1'
	* Penalisation
	local j 1
	matrix `P1'=J(1,`c',0)
	while `j'<`c' {
        scalar `bj'=`b'[1,`j']
		scalar `s'=sd[1,`j']
		if (`bj'*`s')>$delta {
			matrix `P1'[1,`j']=$lambda*`s'
		}
		else if (`bj'*`s')<-$delta {
			matrix `P1'[1,`j']=-$lambda*`s'
		}
		else {
			matrix `P1'[1,`j']=$lambda*((3*`s'^2*`bj')/(2*$delta)-(`s'^4*`bj'^3)/(2*$delta^3))
		}
		local j=`j'+1
	}
	matrix `g' = `g' - `P1'
	if `todo'==1 | `lnf'==. exit

	tempvar d2
	* Second derivatives
	quietly gen double `d2'=exp(`theta')/(1+exp(`theta'))^2
	mlmatsum `lnf' `negH' = `d2'
	*matrix `i' = `negH'
	* !! PR bug fix
	matrix negH = `negH'
	* 2nd derivs of penalisation function zero except in fudge region
	local j 1
	matrix `P2'=J(`c',`c',0)
	while `j'<`c' {
        scalar `bj'=abs(`b'[1,`j'])
		scalar `s'=sd[1,`j']
		if (`bj'*`s')<$delta {
			matrix `P2'[`j',`j']=$lambda*((3*`s'^2)/(2*$delta)-(3*`s'^4*`bj'^2)/(2*$delta^3))
		}
		local j=`j'+1
	}
	matrix `negH' = `negH' + `P2'
end
