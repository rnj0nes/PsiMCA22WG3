* Survival model with covariates
cap confirm file model5mi.out
if _rc~=0 {
   cap erase inpfiles.dat 
   rdoc init inpfiles.dat 
   * run individually but only for purposes of saving inputdatafile 
   forvalues i=1/5 {
      cap erase inp`i'.dat 
      runmplus tc py a1cspline `covrange' if _Imputation_==`i' , ///
      	variable(SURVIVAL = py; TIMECENSORED = tc (0 = NOT 1 = RIGHT);) ///
         model(py on a1cspline `covrange' ; ) ///
         log(off) savelog(model5`i') saveinp(model5`i') saveinputdatafile(inp`i') //
      r inp`i'.dat //`
   }
   rdoc close 
   cap erase model5mi.inp 
   filefilter model51.inp model5mi.inp , from("FILE = inp1.dat ;") to("FILE = inpfiles.dat ; \n TYPE = imputation ;")
   *** type model5mi.inp 
   !$mplus_path $cwd/model5mi.inp
}
runmplus , po(model5mi.out) log(off)
mat E5 = r(estimate)
mat SE5 = r(se)
mat z5 = r(z)
*** runmplus_show_output_segment model5mi.out "MODEL RESULTS" 120
