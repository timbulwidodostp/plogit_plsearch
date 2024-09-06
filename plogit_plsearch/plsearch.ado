*! v1.0.1 25June01 GA
* Grid search routine for plogit
program define plsearch, eclass
	version 6.0
	gettoken plvlist 0:0, parse(",")
	syntax [, LAMbda(numlist asc min=0) LASso noPLot noSTop TRace NOIsy /*
		*/ PRev INIt(str) *]
	* set up locals to switch between AIC/GCV
	local result=cond(missing("`lasso'"),"aic","gcv")
	local sign=cond(missing("`lasso'"),">","<")
	tempvar touse
	tempname `result' gcv_0 
	if missing("`plot'") {
		tempvar x y
		qui gen `x'=.
		qui gen `y'=.
		lab var `x' "lambda"
		lab var `y' "`result'"
	}
	cap est drop _plogit	
	local trace=cond("`trace'"=="","*","")
	if missing("`noisy'") { local qui quietly }
	if missing("`lambda'") { local lambda 0 }
	if !missing("`init'") {
		if "`init'"=="prev" {
			local prev1 prev
			local init
		}
		else local init  init(`init')
	}
	`trace' di in gr "Lambda" _col(20) "`result'" _col(40) "EDF" _col(60) "s"
	gettoken lam lambda:lambda
	if !`lam' {
		* Fit unpenalised model using plogit
		`qui' plogit `plvlist', lam(0) /* fit penalised model */
	}
	else {
		* Fit first plogit model
		`qui' plogit `plvlist', `options' lam(`lam') `init' `prev1' `lasso'
	}
	`trace' di in gr `lam' _col(20) %6.3f `e(`result')' _col(40) %5.2f `e(edf_m)' /*
		*/ _col(60) %3.2f `e(s)'
	if missing("`plot'") {
		qui replace `x'=`lam' in 1
		qui replace `y'=e(`result') in 1
	}
	scalar ``result''=e(`result')
	local init
	if !missing("`prev'") {
		/* use as initial values for next model */
		tempname eb		
		matrix `eb'=e(b)
		local init init(\`eb')
	}
	* Save results to fit null model
	scalar `gcv_0'=-2*e(ll_0)/(e(N)*(1-1/e(N))^2)
	local depvar `e(depvar)'
	gen byte `touse'=e(sample)
	est hold _plogit	
	* Continue iteration until best model found
	local i 1
	local done 0
	while !`done'&!missing("`lambda'") {
		local i=`i'+1
		gettoken lam lambda:lambda
		`qui' plogit `plvlist', lambda(`lam') `options' `lasso' `init'
		`trace' di in gr `lam' _col(20) %6.3f `e(`result')' _col(40) /*
			*/ %5.2f `e(edf_m)' _col(60) %3.2f `e(s)'
		if missing("`plot'") {
			qui replace `x'=`lam' in `i'
			qui replace `y'=e(`result') in `i'
		}
		if !missing("`prev'") { matrix `eb'=e(b) }
		* Update best model
		if `e(`result')'`sign'``result'' {
			scalar ``result''=`e(`result')'
			est drop _plogit
			est hold _plogit 
		}
		else if missing("`stop'") { local done 1 }
	}
	if missing("`plot'") {
		local ymax=``result''
		gr `y' `x', xlab(`lam') ylab(`ymax')
	}
	* Test if null model is better
	local null 0
	if !`done' {
		tempname res_0
		scalar `res_0'=cond(missing("`lasso'"),0,`gcv_0')
		`trace' di in gr "inf" _col(20) %6.3f `res_0' _col(40) %5.2f 0 /*
			*/ _col(60) %3.2f 0
		if `res_0'`sign'``result'' {
			di in bl "Warning: null model may be optimal"
			`qui' plogit `depvar' if `touse', lam(0)
            estimates scalar lambda = -1
			est drop _plogit
			est hold _plogit
		}
	}
	est unhold _plogit
	plogit
end

