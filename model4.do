* Joint model
cap confirm file model4mi.out
if _rc~=0 {
   cap erase inpfiles.dat 
   rdoc init inpfiles.dat 
   * run individually but only for purposes of saving inputdatafile 
   forvalues i=1/5 {
      cap erase inp`i'.dat 
      runmplus tc py rtr200-rtr205 a1cspline `covrange' if _Imputation_==`i' , ///
      	variable(SURVIVAL = py; TIMECENSORED = tc (0 = NOT 1 = RIGHT);) ///
         model(py on a1cspline `covrange' ; ///
            i s | rtr200@0 rtr201@2 rtr202@4 rtr203@6 rtr204@8 rtr205@10; ///
            i s on a1cspline `covrange'; ///
            py on i s; ) ///
         log(off) savelog(model4`i') saveinp(model4`i') saveinputdatafile(inp`i') //
      r inp`i'.dat //`
   }
   rdoc close 
   cap erase model4mi.inp 
   filefilter model41.inp model4mi.inp , from("FILE = inp1.dat ;") to("FILE = inpfiles.dat ; \n TYPE = imputation ;")
   type model4mi.inp 
   !$mplus_path $cwd/model4mi.inp
}
runmplus , po(model4mi.out) log(off)
mat E4 = r(estimate)
mat SE4 = r(se)
mat z4 = r(z)
runmplus_show_output_segment model4mi.out "MODEL RESULTS" 205
