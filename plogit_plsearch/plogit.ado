*! v 1.0.8 ARB/GA/PR 20jul2007
* Translated to Stata 8.2 PR 20/12/2005
* 1.0.8: noDrop option permanently on to combat conformability problems with lasso
* 1.0.7: Fixed problem with loglikelihood programs being r-class
* 1.0.6: Fixed problem with collinear variables & perefect predictions
* 1.0.4: Adds refining estimates and new initial values options
* 1.0.3: Fixes problem in 1.0.2 with collinear variables causing conformability error
* Helper files (must be installed):
*	plogi_l1.ado
*	plogi_l2.ado (1.0.4: bug in fudge region derivatives fixed)
program define plogit, eclass
	version 8.2
	if replay() {
		if "`e(cmd)'" != "plogit" {
			error 301
		}
		Replay `0'
	}
	else {
		Estimate `0'
	}
end

program define Estimate, eclass
	version 8.2
	syntax varlist(min=1) [if] [in] [fweight pweight] [, /*
		*/ LAMbda(real 0) noLOg or Level(passthru) CHange LASso /*
		*/ DELta(real 1e-5) EPSilon(real 1e-5) PREV INIt(passthru) /*
		*/ noREFine noDROP ]
	marksample touse
	markout `touse', strok
	gettoken yvar varlist:varlist
	qui inspect `yvar' if `touse'
	if r(N_unique)==1 {
		di as err "Outcome does not vary"
		exit 2000
	}
	else if r(N_unique)!=2 {
		di as err "Outcome not binary"
		exit 2000
	}
	global lambda `lambda'
	global delta `delta'
	local drop nodrop // GA: bug fix for lasso
	if !missing("`log'") {
		local qui quietly
	}
	* Starting values 1
	if !missing("`prev'") {
		* Take starting values from previous model
		if !missing("`init'") {
			di in bl "Warning: prev option ignored"
		}
		else {
			tempname beta
			matrix `beta'=e(b)
			local init init(`beta',skip)
		}
	}
	* Drop collinear variables and report
	noi _rmcoll `varlist'
	local expltmp `r(varlist)' /* temp varlist */
	* Store unpenalised ML estimates and vector length
	qui logit `yvar' `expltmp' [`weight'`exp'] if `touse'
	tempname obeta olength openal oprob
	matrix `obeta'=e(b)
	local vl = colsof(`obeta')-1	/* Note excludes constant (always last element in beta vector) */
	* Detect whether any vars perfectly predict response
	qui count if `touse'
	if `r(N)'>`e(N)' {
		tempname names
		if `vl' {
			matrix `names' = `obeta'[1,1..`vl']
			local explvar : colnames(`names')
		}
		tokenize `expltmp'
		while !missing("`1'") {
			local found 0
			local i 0
			while `i'<`vl'&!`found' {
				local i=`i'+1
				local explvri: word `i' of `explvar'
				if "`1'"=="`explvri'" {
					local found 1
				}
			}
			if !`found' { 
				noi di as txt "`1' predicts outcome perfectly: dropped"
			}
			mac shift
		}
		local ndrop=r(N)-e(N)
		if `ndrop'>1 local s s
		noi di as txt "`ndrop' observation`s' dropped"
		qui replace `touse'=e(sample)
	}
	else local explvar `expltmp'
	* Penalisation matrix and SD vector (not including the constant):
	matrix P=I(`vl'+1)
	matrix mean=vecdiag(P)
	matrix sd=vecdiag(P)
	tokenize `explvar'
	local i 0
	while `i'<`vl' {
		local i=`i'+1
		qui summ ``i''
		matrix P[`i',`i']=r(Var)
		matrix mean[1,`i']=r(mean)
		matrix sd[1,`i']=r(sd)
	}
	matrix P[`vl'+1,`vl'+1]=0
	matrix sd[1,`vl'+1]=0
	* Calculate original length of unpenalised beta vector exc. _cons (needed for scaled lasso parameter)
	local i 1
	scalar `olength'=0
	scalar `openal'=0
	while `i' <= `vl' {
		scalar `olength'=`olength'+abs(`obeta'[1,`i']*sd[1,`i'])
		scalar `openal'=`openal'+(`obeta'[1,`i']*sd[1,`i'])^2
		local i=`i'+1
	}
	* Starting values 2
	if missing("`prev'")&missing("`init'") {
		* Create penalised starting values from either null or unpenalised model
		if missing("`lasso'") {
			// Choose either null or unpenalised based on penalised ll
			if ($lambda*`openal'<e(chi2)) {
				local init init(`obeta')
			}
		}
		* Create lasso starting values from penalised model
		else {
			global lambda=2*`lambda'*`olength'/`openal' /* Equivalent lambda */
			// Create starting values for penalised model
			if ($lambda*`openal'>e(chi2)) {
				// Start with null values instead of unpenalised values
				qui logit `yvar' [`weight'`exp'] if `touse'
			}
			`qui' di as txt "Obtaining starting values:" _c
			`qui' ml model d2 plogi_l1 (`yvar'=`explvar') [`weight'`exp'] if `touse', /*
				*/ missing continue maximize search(off) nopreserve
			global lambda=`lambda'
			tempname beta
			matrix `beta'=e(b)
			local b0 `beta'[1,`vl'+1]
			local i 0
			while `i'<`vl' {
				* Seems to provide better starting values?
				local i=`i'+1
				if abs(`beta'[1,`i']*sd[1,`i'])<0.05 {
					matrix `b0'=`b0'+`beta'[1,`i']*mean[1,`i']
					matrix `beta'[1,`i']=0
				}
			}
			local init init(`beta')
		}
	}
	if missing("`lasso'") {
		* Penalised
		* Calculate loglik for constant-only model - needed for model LR
		qui logit `yvar' [`weight'`exp'] if `touse'
		`qui' ml model d2 plogi_l1 (`yvar'=`explvar') [`weight'`exp'] if `touse', /*
			*/ missing continue maximize `init' search(off) nopreserve /*
			*/ title(Penalised logistic regression)
		tempname edf IV
		matrix `IV'=negH*e(V)
		*matrix `IV'=r(I)*e(V)
		scalar `edf'=trace(`IV')-1	/* edf for penalised model */
	}
	else {
		* Lasso
		`qui' di as txt _n "Estimating model:" _c
		qui logit `yvar' [`weight'`exp'] if `touse'
		`qui' ml model d2 plogi_l2 (`yvar'=`explvar') [`weight'`exp'] if `touse', /*
			 */ missing continue maximize `init' search(off) nopreserve /*
			 */ title("Penalised logistic regression (lasso)")
		tempname edf newin
		if missing("`refine'") {
			local its 0
			local done 0
			while !`done' {
				local its=`its'+1
				`qui' di as txt _n "Refining estimates(`its'):" _c
				matrix `newin'=e(b)
				local b0 `newin'[1,`vl'+1]
				local comma
				local j 1
				local i 1
				while `i'<=`vl' {
					gettoken vari explvar:explvar
					if abs(_b[`vari']*sd[1,`i'])>=`epsilon' {
						local explvar `explvar' `vari'
						if missing("`drop'") {
							matrix mean[1,`j']=mean[1,`i']
							matrix sd[1,`j']=sd[1,`i']
							local j=`j'+1
						}
					}
					else {
						matrix `b0'=`b0'+_b[`vari']*mean[1,`i']
						if missing("`drop'") {
							if missing("`comma'") { local newline _n }
							`qui' di as txt `newline' "`comma'`vari'" _c
							local comma ", "
							local newline
						}
						else {
							local explvar `explvar' `vari'
							matrix `newin'[1,`i']=0
						}
					}
					local i=`i'+1
				}
				if !missing("`comma'") {
					`qui' di as txt " dropped"
				}
				if missing("`drop'") {
					matrix mean=mean[1,1..`j']
					matrix mean[1,`j']=1
					matrix sd=sd[1,1..`j']
					matrix sd[1,`j']=0
					local vl=`j'-1
				}
				* Calculate loglik for constant-only model - needed for model LR
				qui logit `yvar' [`weight'`exp'] if `touse'
				`qui' ml model d2 plogi_l2 (`yvar'=`explvar') [`weight'`exp'] if `touse', /*
					 */ missing continue maximize init(`newin',skip) search(off) nopreserve /*
					 */ title("Penalised logistic regression (lasso)")
				* Check if any betas<epsilon: only do once if nodrop specified
				local done 1
				local i 0
				while `i'<`vl'&`done'&missing("`drop'") {
					local i=`i'+1
					local vari: word `i' of `explvar'
					if abs(_b[`vari']*sd[1,`i'])<`epsilon' {
						local done 0
					}
				}
 			}
		}
	}
	ereturn scalar ll_unp = $PL_llnp
	ereturn scalar ll_penal = $PL_llpen
	ereturn scalar chi2_unp = -2*(e(ll_0)-$PL_llnp)
	* Calculate edf for lasso model & stick in `edf'
	if !missing("`lasso'") edf_las `edf' `epsilon'
	ereturn scalar edf_m = `edf'
	ereturn local vars `varlist'
	ereturn scalar lambda = $lambda
	ereturn scalar aic = e(chi2_unp)-2*e(edf_m)
	ereturn scalar gcv = -2*e(ll_unp)/(e(N)*(1-(e(edf_m)+1)/e(N))^2)
	ereturn local cmd "plogit"
	tempname pbeta plength
	matrix `pbeta'=e(b)
	local vl=colsof(`pbeta')
	local i 1
	scalar `plength'=0
	while `i' < `vl' {		/* Note excludes constant (always last element in beta vector) */
		scalar `plength'=`plength'+abs(`pbeta'[1,`i']*sd[1,`i'])
		local i=`i'+1
	}
	ereturn scalar s = `plength'/`olength'
	Replay, `level' `or'
	if !missing("`change'") {
		tempname pbetai
		di as txt _col(10) "|  Change (compared to ML estimate)"
		local vn: colnames `obeta'
		tokenize `vn'
		local i 1
		while !missing("`1'") {
			if "`1'" != "_cons" {
				cap scalar `pbetai'=_b[`1']
				if _rc {
					scalar `pbetai'=0
				}
				di as txt %8s "`1'" " | " as res %8.1f 100*(`pbetai'/`obeta'[1,`i']-1) "%"
			}
			local i=`i'+1
			mac shift
		}
		di as txt "{hline 78}"
	}
end

program define Replay
	version 8.2
	syntax [, Level(int $S_level) or noml]
	if !missing("`or'") {
		local or eform(Odds ratio)
	}
	if missing("`ml'") {
		di _n as txt "`e(title)'" _col(51) "Number of obs   =" as res %11.0f `e(N)'
		di as txt _col(51) "Effective df    =" as res %11.2f `e(edf_m)'
		di as txt "log likelihood = " as res %10.0g `e(ll)' /*
			*/ _col(51) as txt "LR chi2 (no pen)=" as res %11.2f `e(chi2_unp)'
		ml display, level(`level') `or' `plus' noheader
	}
end

*! GA 8May01
* Calculates edf of lasso fit
program define edf_las
	version 8.2
	args edf_las epsilon
	tempvar xb p d dxixj trace
	tempname coefs cons lambdaW XtDX I Iinv XIinvXD edf mean sd
	* obtain info from plogit fit+process options
	qui predict `xb' if e(sample)
	qui gen `p'=1/(1+exp(-`xb'))
	qui gen `d'=`p'*(1-`p')
	* obtain columns of X: X_j & calculate lambdaW
	gen byte `cons'=1
	matrix `coefs'=e(b)
	local vars:colnames(`coefs')
	local nvars : word count `vars'
	matrix `mean'=mean
	matrix `sd'=sd
	local nj 1
	local j 1
	while `j'<`nvars' {
		gettoken varj vars:vars
		if abs(_b[`varj']*sd[1,`j'])>=`epsilon' {
			tempname X_`nj'
			qui gen `X_`nj''=(`varj'-mean[1,`j'])/(sd[1,`j'])
			local var`nj' `varj'
			matrix `mean'[1,`nj']=mean[1,`j']
			matrix `sd'[1,`nj']=sd[1,`j']
			local nj=`nj'+1
		}
		local j=`j'+1
	}
	local X_`nj' `cons'
	matrix `sd'=`sd'[1,1..`nj']
	matrix `sd'[1,`nj']=0
	matrix `lambdaW'=J(`nj',`nj',0)
	local j 1
	while `j'<`nj' {
		matrix `lambdaW'[`j',`j']=$lambda/(abs(_b[`var`j'']*`sd'[1,`j']))
		local j=`j'+1
	}
	* calculate inv(X'DX+lambdaW)
	matrix `lambdaW'[`nj',`nj']=0
	gen `dxixj'=0 
	local i 1
	matrix `XtDX'=J(`nj',`nj',0)
	local i 0
	while `i'<`nj' {
		local i=`i'+1
		local j 0
		while `j'<`nj' {
			local j=`j'+1
			qui replace `dxixj'=`d'*`X_`i''*`X_`j''
			qui summ `dxixj', meanonly
			matrix `XtDX'[`i',`j']=r(sum)
		}
	}
	matrix `I'=`XtDX'+`lambdaW'
	matrix `Iinv'=syminv(`I')
	* calculate trace elements
	scalar `edf'=0
	qui gen `XIinvXD'=.
	local j 0
	while `j'<`nj' {
		local j=`j'+1
		qui replace `XIinvXD'=`X_`j''*`X_`j''*`Iinv'[`j',`j']*`d'
		summ `XIinvXD', meanonly
		scalar `edf'=`edf'+r(sum)
		local k `j'
		while `k'<`nj' {
			local k=`k'+1
			* Calculate edf for lasso model & stick in `edf'
			qui replace `XIinvXD'=`X_`j''*`X_`k''*`Iinv'[`j',`k']*`d'
			summ `XIinvXD', meanonly
			scalar `edf'=`edf'+2*r(sum)
		}
	}
	scalar `edf_las'=`edf'-1	/* remove constant */
end
