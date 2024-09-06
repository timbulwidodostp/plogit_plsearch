*! GA v1.1.6 25/07/03
program define stepdown, rclass
	version 7.0
	local qui quietly
	if missing("`e(cmd)") { error 301 }
	* DISENTANGLE CLUSTERS - most of this adapted from -sw-
	local nclust 0 /* # clusters */
	gettoken tok : 0, parse(" ,[")
	IfEndTrm "`tok'"
	while `s(IsEndTrm)'==0 {
		gettoken tok 0 : 0, parse(" ,[")
		if substr("`tok'",1,1)=="(" {
			local list
			while substr("`tok'",-1,1)!=")" {
				if "`tok'"=="" { 
					di in red "varlist invalid"
					exit 198
				}
				local list "`list' `tok'"
				gettoken tok 0 : 0, parse(" ,[")
			}
			local list "`list' `tok'"
			unabbrev `list'
			local nclust = `nclust' + 1
			local clust`nclust' `s(varlist)' /* vars in cluster */
			local n`nclust' `s(k)' /* # vars in cluster */
			local vars `vars' `s(varlist)'
		}
		else {
			unabbrev `tok'
			local i 1
			local w : word 1 of `s(varlist)'
			while !missing("`w'") {
				local nclust=`nclust' + 1
				local clust`nclust' `w'
				local n`nclust' 1
				local i = `i' + 1
				local w : word `i' of `s(varlist)'
			}
			local vars `vars' `s(varlist)'
		}
		gettoken tok : 0, parse(" ,[")
		IfEndTrm "`tok'"
	}
	local nvars: word count `vars'
	* SYNTAX
	syntax [if] [in] [,r2(real 0) LOCKterm1 noFIll TRace HOld SLOW]
	marksample touse
	local trace=cond(missing("`trace'"),"*","")
	if `r2'>1|`r2'<0 {
		noi di in red "r2 should lie between 0 and 1"
		exit 198
	}
	if missing("`fill'") { /* fill in "gaps" in varlist */
		tempname coefs
		matrix `coefs'=e(b)
		local mvars: colnames(`coefs') /* model variables */
		local nmvars `e(df_m)'
		local i 0
		while `i'<`nmvars' {
			local i=`i'+1
			gettoken mvar`i' mvars:mvars
			local in`i' 0 /* indicates if variable is in `vars' */
		}
		local i 0
		while `i'<`nvars' { /* see if `vars' is in `mvars' */
			local i=`i'+1
			local vari: word `i' of `vars'
			local done 0
			local j 0
			while `j'<`nmvars' {
				local j=`j'+1
				if "`vari'"=="`mvar`j''" {
					local in`j' 1
					local done 1
				}
			}
			if !`done' {
				di in red "`vari' not in model"
				exit 198
			}
		}
		local i 0
		while `i'<`nmvars' { 
			local i=`i'+1
			if !`in`i'' {
				local nclust=`nclust'+1
				local clust`nclust' `mvar`i''
				local n`nclust' 1
				local nvars=`nvars'+1
				local vars `vars' `mvar`i''
			}
		}
		if `nvars'>`nmvars' {
			di in red "repeated predictors in varlist"
			exit 198
		}
	}
	* SETUP
	tempvar xb
	tempname r2_cur r2_min r2_new r2_tmp tss rss df_r
	`qui' predict `xb', xb
	local weight `e(wtype)'
	local exp \`e(wexp)'
	if !missing("`hold'") {
		cap estimates drop _stepdwn
		if !_rc { di in bl "Warning: previous estimates found & deleted" }
		cap est hold _stepdwn, restore
		if _rc {
			est hold _stepdwn /* version 6 does not have restore option */
			local v6hold 1
		}
	}
	`qui' regress `xb' `vars' if `touse' [`weight'`exp']
	scalar `rss'=e(rss)
	scalar `df_r'=e(df_r)
	scalar `tss'=e(mss)+e(rss)
	if `e(df_m)'<`nvars' {
		di in bl "Warning: predictors have been dropped - collinearity?" _n
	}
	scalar `r2_cur'=e(r2)
	local i 0
	while `i'<`nclust' {
		local i=`i'+1
		local flag`i' 1 /* indicates whether to consider cluster */
		local in`i' 1 /* indicates whether in model */
	}
	if !missing("`lockterm1'") { local flag1 0 }
	local done 0
	local its 0
	if `r2_cur'<`r2' {
		di in bl "Warning: R-squared is already lower than `r2'"
		local done 1
	}
	if abs(`r2_cur'-1)>5e-4 {  /* r2_cur <=.0.999 */
		di in bl "Warning: initial R-squared value less than 1" _n
	}
	else di in gr _col(23) "begin with full model" _col(70) "R2 = " in ye /* 
	 */ %4.3f `r2_cur'
	* MAIN LOOP
	while !`done' {
		local its=`its'+1
		scalar `r2_min'=-1
		local besti 0
		`trace' di
		local i 0
		while `i'<`nclust' {
			local i=`i'+1
			if `in`i'' & `flag`i'' {	/* loop over clusters left in model */
				if missing("`trace'") {
					local drop `clust`i''
					if length("`drop'")>45 { 
						local drop=substr("`drop'",1,42)
						local drop `drop'...
					}
				}
				* use slow option whenever R2 is v.close to 1
				if abs(`r2_cur'-1)>5e-4&missing("`slow'") {
					`qui' test `clust`i''	
					scalar `r2_new'=1-`rss'/`tss'*(1+r(F)*`n`i''/`df_r')
				}
				else {
					local var_cur
					local j 0
					while `j'<`nclust' {
						local j=`j'+1
						if `j'!=`i' & `in`j'' { 
							local var_cur `var_cur' `clust`j''
						}
					}
					`qui' regress `xb' `var_cur' if `touse' [`weight'`exp']
					scalar `r2_new'=e(r2)
				}
				`trace' scalar `r2_tmp'=round(`r2_new',1e-8) /* to avoid -0.000 */
				`trace' di in gr "drop in R2 = " in gr  %4.3f (`r2_cur'-`r2_new') /*
				  */ in gr _col(23) "`drop'" _col(70) "R2 = " in gr /*
				  */%4.3f `r2_tmp' 
				if (`r2'-`r2_new')>1e-8 { local flag`i' 0 }	   /* var too important to drop */
				else {
					if `r2_new'>`r2_min' {
						scalar `r2_min'=`r2_new'
						local besti `i'
					}
				}	
			}
		}
		* DISPLAY
		if `besti' {
			local in`besti' 0
			local drop `clust`besti''
			if length("`drop'")>36 {
				local drop=substr("`drop'",1,33)
				local drop `drop'...
			}
			if missing("`slow'") { /* provide estimates for next iteration */
				local var_cur
				local i 0
				while `i'<`nclust' {
					local i=`i'+1
					if `in`i'' { local var_cur `var_cur' `clust`i'' }
				}
				`qui' regress `xb' `var_cur' if `touse' [`weight'`exp']
				scalar `rss'=e(rss)
				scalar `df_r'=e(df_r)
				scalar `r2_min'=e(r2) /* replace analytic value */
			}
			di in gr "drop in R2 = " in ye %4.3f (`r2_cur'-`r2_min') _col(23) /*
			 */ in gr "removing `drop'" _col(70) "R2 = " in ye %4.3f /*
			 */ `r2_min'
			scalar `r2_cur'=`r2_min'
		}
		else local done 1
	}
	* FINAL FIT
	local vars
	local nclust2 0
	local i 0
	while `i'<`nclust' {
		local i=`i'+1
		if `in`i'' { 
			local nclust2=`nclust2'+1
			local vars `vars' `clust`i''
		}
	}
	local nvars: word count `vars' 
	regress `xb' `vars' if `touse' [`weight'`exp'], depname(mu) 
	if !missing("`v6hold'") { est unhold _stepdwn }
	ret scalar nclust=`nclust2' /* # clusters in final model */
	ret scalar nvars=`nvars' /* # vars in final model */
	ret local vars `vars' /* final varlist */
end


program define IfEndTrm, sclass /* taken from -sw- */
	sret local IsEndTrm 1
	if "`1'"=="," | "`1'"=="in" | "`1'"=="if" | /*
	*/ "`1'"=="" | "`1'"=="[" { exit }
	sret local IsEndTrm 0
end
