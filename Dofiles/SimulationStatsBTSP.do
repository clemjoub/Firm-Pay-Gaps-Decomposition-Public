clear
clear matrix
clear mata
set more off
set trace off
set matsize 4000
timer clear 1
timer on 1
set seed 1

********************************************************************************
* Firm Pay gaps : a decomposition
*
*	This do file estimates the parameters used as input into the fortran 
*	simulation code
*
*	Input: sample.dta (data format is one observation per employee and year)
* 	
*	Variables used:
*
* 	upi - unique employee identifier
* 	year - year in which salary is observed
* 	active - currently employed
* 	salarycurr - annual nominal salary
* 	lagsal - annual nominal salary in previous year
* 	logincr - log current salary - log of last year salary
* 	noincrease - dummy for no change in nominal salary
* 	CPI - CPI index
* 	newGspell - dummy for first year in which employee is observed at current 
*				grade (we ignore ungraded employment spells)
* 	endGspell - dummy for last year in which employee is observed at current grade
* 	Gtenure - number of years spent in current spell
* 	firstGyear - first year in which employee was observed in a graded position
* 	type - category of employee (gender*grade* nationality group)

*	Notes on timing:
*	
*	salaries are considered 'end of period'. Individuals in their last year are 
*	excluded from the 'end of period' distribution and included in the 'attriter'
*	distribution. The model simulations start in 1988, using 1987 'end of period'
*	salaries as the distribution of salaries before entries,exits and wage growth.
*
********************************************************************************

*set parameters
	global nbboot = 1				// if bootstrapping input parameter estimates		
	global ktp=24					// total number of groups of employees (entry grade * gender * nationality group)
	global maxsal=500.0				// maximum salary considered
	global binsize=1.0				// size of salary bins	
	global nbsalbins=$maxsal/$binsize
	global ct1=1987					// first year of data
	global ctT=2015					// last year of data
	
* set paths
	global path="" //set project path
	global datapath="" //set path to personnel records
	global outpath="$path\Estimation\AggregateGap"

use "$datapath\sample.dta", clear

	keep upi year active salarycurr lagsal  logincr noincrease CPI  newGspell endGspell Gtenure firstGyear type  
	keep if active ==1 

	
*Set exit indicator to missing in last year
*	replace endspell=. if year==${ctT}
	
*Some variable definitions
	quietly do "$path/Dofiles/Simstat/VariableDefs.do"
	save "$path\temp", replace

	
*loop over bootstrap replications
	forvalues boot=1(1)$nbboot {
		global j=0
		
		use "$path\temp", clear 
		if `boot'>1 {   //first set of estimates is without bootstrapping
			bsample		
			display "bootstrap # `boot'"
		} //if `boot'
		
		
		forvalues tp=1(1)$ktp{ //loop over groups of employees
			forvalues i=$ct1(1)$ctT{  //loop over years of data
					preserve	
						global condition1="ct==`i' & type==`tp' " //keep records corresponding to a group and a year
						keep if $condition1
						global yr=`i'
						global j=$j+1					
						
						quietly do "$path/Dofiles/Simstat/SimstatEntry.do" 		//entry parameters
						quietly do "$path/Dofiles/Simstat/SimstatNoincrease.do" //salary growth parameters 1			
						quietly do "$path/Dofiles/Simstat/SimstatGrowth.do"		//salary growth parameters 2
						quietly do "$path/Dofiles/Simstat/SimstatExits.do"		//exit parameters 1
						quietly do "$path/Dofiles/Simstat/SimstatEndOfPeriod.do"	//exit parameters 2

					restore
			} //forvalues i
		} //forvalues tp
			
		*Store results in a matrix	
		mat simstats=(entry, noincrease,growth,exit,eop, out, inter)
		mat colnames simstats =  nentry muentry sdentry noincrease  betagrowth betagrowth_se alphagrowth siggrowth betaexit beta2exit alphaexit nend muend sigend nout muout sigout nint muint sigint
		mat2txt, m(simstats) sav("$outpath\simstats`boot'") format(%15.9f) replace
		xsvmat simstats, saving("$outpath\simstats.dta", replace) names(col)
		
		timer off 1
		timer list 1
	} //forvalues boot

	*Produce salary histograms for number of entrants, leavers and for the end of period stock
	preserve
		keep if ct<$ctT+1 & ct>$ct1-1
		collapse (sum) exits (sum) entries (sum) eop if active==1, by(type ct salbin)
		sort salbin type ct
		outfile salbin type ct entries exits  eop using "$outpath\salhist.txt" , nol replace wide
	restore


clear
