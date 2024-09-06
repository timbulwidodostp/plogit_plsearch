*! v 1.0.1 ARB/PR  21dec2005
program define plogi_l1
	version 8.2
	args todo b lnf g negH
	tempvar theta lnfj
	tempname bPb penal P1 P2 llnp
	mleval `theta'=`b'
	quietly gen double `lnfj'=-ln(1+exp(-`theta')) if $ML_y1==1
	quietly replace `lnfj'=-ln(1+exp(`theta')) if $ML_y1==0
	mlsum `lnf'=`lnfj'
	scalar `llnp' = `lnf'
	* Penalisation
	matrix `bPb'=`b'*P*`b''
	scalar `penal'=0.5*$lambda*`bPb'[1,1]
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
	matrix `P1'=($lambda*P*`b'')'
	matrix `g' = `g' - `P1'
	if `todo'==1 | `lnf'==. exit

	tempvar d2
	* Second derivatives
	quietly gen double `d2'=exp(`theta')/(1+exp(`theta'))^2
	mlmatsum `lnf' `negH' = `d2'
	* Penalisation
	matrix `P2'=$lambda*P
	* !! PR bug fix
	*matrix `i' = `negH'
	matrix negH = `negH'
	matrix `negH' = `negH' + `P2'
end
