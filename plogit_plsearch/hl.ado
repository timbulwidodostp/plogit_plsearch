*! v3.0.5 23Jan04 GA & ARB
/*
	Calculate HL and c-index
*/
program define hl, rclass sortpreserve byable(recall)
	version 8.0
	syntax [varlist(default=none)] [if] [in] /*
	 */ [, Group(int 10) CUTpoints(numlist) Q(str) CR /*
	 */    SCale(str) RANge ZEro(str) OUTsample /*
	 */    noROC noTABle      /*
	 */    PLot BYPLot(str) FW G1opt(str) G2opt(str) ]

	quietly { /* SOQ */

	* Identify main variables
	if missing("`varlist'") {
		if "`e(cmd)'"=="logit"|"`e(cmd)'"=="logistic"|"`e(cmd)'"=="xtlogit"| /*
		 */ ("`e(cmd)'"=="xtgee"&"`e(family)'"=="binomial"&"`e(link)'"=="logit") {
			tempvar xb prob
			local y `e(depvar)'
			predict `xb', xb
			gen `prob' = 1/(1+exp(-`xb'))
		}
		else {
			error 102
		}
	}
	else {
		tokenize `varlist'
		if missing("`2'") {
			error 102
		}
		if !missing("`3'") {
			error 103
		}
		args y prob
	}

	* Identify estimation sample & check vars
	marksample touse
	markout `touse' `y' `prob'
	capture assert `y'==0 | `y'==1 if `touse'
	if _rc {
		noi di "{error}`y' should be 0/1"
		exit 198
	}

	* option: zero
	if !missing("`zero'") {
		count if !`prob' & `touse'
		if r(N)>0 {
			tempvar prob0
			gen `prob0' = `prob' if `touse'
			recode `prob0' 0 = `zero'
			local prob `prob0'
		}
	}

	capture assert `prob'>=0 & `prob'<=1 if `touse'
	if _rc {
		noi di "{error}`prob' contains values outside of [0,1]"
		exit 198
	}

	* Deal with options
	* option: scale
	if !missing("`scale'")&"`scale'"!="pe"&"`scale'"!="pr"&"`scale'"!="lo" {
		noi di in red "scale(`scale') is invalid"
		exit 198
	}
	* options: priority - q, cr, cutpoints, group
	if missing("`q'") {
		if !missing("`cr'") {
			local cutpoints "0.025 0.05 0.1 0.2 "
		}
		local qgen 1
		cap label drop _hlqlab
		label define _hlqlab -1 "Total"
		tempvar q
		if !missing("`cutpoints'") {
			gen byte `q' = 1 if `touse'
	   		if "`scale'"=="lo" {
				label define _hlqlab 1 "-inf -", add
			}
			else if "`scale'"=="pr" {
				label define _hlqlab 1 ".0 -", add
			}
			else {
				label define _hlqlab 1 "0% -", add
			}
			local i 2
			foreach cut of local cutpoints {
				if `cut'>0 & `cut'<1 {
					replace `q' = `q'+(`prob'>=`cut') if `touse'
			   		if "`scale'"=="lo" {
						local cut = logit(`cut')
						local cut : display %5.3f `cut'
						label define _hlqlab `i' "`cut' -", add
					}
					else if "`scale'"=="pr" {
						local cut : display %5.3f `cut'
						label define _hlqlab `i' "`cut' -", add
					}
					else {
						local cut=`cut'*100
						local cut : display %5.1f `cut'
						label define _hlqlab `i' "`cut'% -", add
					}
					local ++i
				}
			}
		}
		else xtile `q' = `prob' if `touse', n(`group')
		label values `q' _hlqlab
	}
	* Plotting options
	if !missing("`g1opt'")|!missing("`g2opt'")|!missing("`byplot'")|!missing("`fw'") {
		local plot plot
		if !missing("`byplot'") {
			local byopt by(`byplot')
		}
		if !missing("`fw'") {
			local fw "[fw=ni]"
		}
	}

	* Main routine
	preserve
	keep if `touse'
	collapse (sum) oi=`y' ei=`prob' (count) ni=`y' (min) mini=`prob' /*
	 */ (max) maxi=`prob', by(`q' `byplot')
	gen hli = (oi-ei)^2/(ei*(1-ei/ni))

	* Totals
	local df = _N
	local nobs = _N+1
	set obs `nobs'
	foreach var of varlist oi ei ni hli {
		summ `var'
		local `var'sum = r(sum)
		replace `var' = r(sum) in `nobs'
	}
	
	* Scalings - log-odds, proportions, %s
	if !missing("`plot'")|("`table'"!="notable") {
		if !missing("`range'") {
			replace mini = mini[1] in `nobs'
			replace maxi = maxi[`df'] in `nobs'
		}
		if "`scale'"=="lo" {
			gen oi2 = log(oi/(ni-oi))
			gen ei2 = log(ei/(ni-ei))
			if !missing("`range'") {
				replace mini = logit(mini)
				replace maxi = logit(maxi)
			}
			local scname "(lo)"
			local scplot "(log-odds)"
		}
		else {
			local scopt = cond("`scale'"=="pr",1,100)
			gen oi2 = `scopt'*oi/ni
			gen ei2 = `scopt'*ei/ni
			if !missing("`range'")&"`scale'"!="pr" {
				replace mini = 100*mini
				replace maxi = 100*maxi
			}
			local scname = cond("`scale'"=="pr","(pr)","(%)")
			local scplot = cond("`scale'"=="pr","(proportion)","(%)")
		}
	}
	
	* Plot
	if !missing("`plot'") {
		lab var oi2 "observed"
		lab var ei2 "predicted"
		scatter oi2 ei2 `fw' in 1/`df' , `g1opt' || /*
		 */ line ei2 ei2 in 1/`df' , sort `byopt' ytitle("Risk `scplot'") /*
		 */ xtitle("Predicted risk `scplot'") `g2opt'
		 * title("Observed and predicted risk (Hosmer-Lemeshow test)")
	}

	* Table
	if "`table'"!="notable" {
		if missing("`range'") {
			local list "oi oi2 ei ei2"
		}
		else {
			local rname1 Obs
			local rname2 Exp
			local list "oi2 ei2 mini maxi"
		}
		char `q'[varname] "Group"
		char oi[varname] "Obs"
		char ei[varname] "Exp"
		char ni[varname] "N"
		char hli[varname] "HL"
		char mini[varname] "min"
		char maxi[varname] "max"
		char oi2[varname] "`rname1'`scname'"
		char ei2[varname] "`rname2'`scname'"
		format hli %9.1f
		format ei %9.1f
		if "`scale'"=="pr"|"`scale'"=="lo" {
			local fmt %9.3f
		}
		else local fmt %9.1f
		format oi2 `fmt'
		format ei2 `fmt'
		format mini `fmt'
		format maxi `fmt'
		if !missing("`qgen'") {
			* q was created
			recode `q' .=-1
		}
		else {
			* q was given - could be str or numeric
			cap confirm string variable `q'
			if _rc {
				* numeric
				tempvar qs
				cap decode `q', gen(`qs')
				if _rc {
					* q didn't have value labels
					gen `qs'=string(`q')
				}
				local q `qs'
				char `qs'[varname] "Group"
			}
			replace `q'="Total" in `nobs'
		}
		noi list `q' `byplot' ni `list' hli , noobs subvar sep(`df')
	}

	} /* EOQ */

	* Output
	di _n _col(8) "{text}number of observations = {result}" %9.0g `nisum'
	di _col(14) "{text}number of groups = {result}" %9.0g `df'
	if missing("`outsample'") {
		local df = `df'-2
	}
	tempname pvalue
	scalar `pvalue' = chiprob(`df',`hlisum')
	di _col(7) "{text}Hosmer-Lemeshow chi2({result}" `df' /*
	 */ "{text}) = {result}" %12.2f `hlisum' _n /*
	 */ _col(19) "{text}Prob > chi2 = {result}" %14.4f `pvalue'

	restore

	* ROC area
	if "`roc'"!="noroc" {
		tempname min max x nn mm ROC
		qui ranksum `prob' if `touse', by(`y')
		local V1 = r(group1)
		local W = r(sum_obs)
		qui summ `y' if `touse'
		scalar `min' = r(min)
		scalar `max' = r(max)
		if `V1' == `min' {
			scalar `x' = `min'
			scalar `min' = `max'
			scalar `max' = `x'
		}
		qui summ `prob' if `touse' & `y'==`max'
		scalar `nn' = r(N)
		qui summ `prob' if `touse' & `y'==`min'
		scalar `mm' = r(N)
		scalar `ROC' = (`nn'*`mm' + `nn'*(`nn'+1)/2 - `W')/(`nn'*`mm')
		di _col(22) "{text}ROC area = {result}" %14.4f `ROC'
		ret scalar roc = `ROC'
	}

	ret scalar hl = `hlisum'
	ret scalar df = `df'
	ret scalar p = `pvalue'
end

