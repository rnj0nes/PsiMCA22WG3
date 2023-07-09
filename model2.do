cap confirm file model2mi.out 
if _rc~=0 {
   cap erase inpfiles.dat 
   rdoc init inpfiles.dat 
   * run individually but only for purposes of saving inputdatafile 
   forvalues i=1/5 {
      cap erase inp`i'.dat 
      runmplus rtr200-rtr205 a1cspline if _Imputation_==`i' , est(MLR) ///
         model(i s | rtr200@0 rtr201@2 rtr202@4 rtr203@6 rtr204@8 rtr205@10; ///
         i s on a1cspline ;) ///
         log(off) savelog(model2`i') saveinp(model2`i') saveinputdatafile(inp`i') //`
      r inp`i'.dat //`
   }
   rdoc close 
   cap erase model2mi.inp 
   filefilter model21.inp model2mi.inp , from("FILE = inp1.dat ;") to("FILE = inpfiles.dat ; \n TYPE = imputation ;")
   type model2mi.inp 
   !$mplus_path $cwd/model2mi.inp
}
runmplus , po(model2mi.out) log(off)
mat E2 = r(estimate)
mat SE2 = r(se)
mat z2 = r(z)
runmplus_show_output_segment model2mi.out "MODEL RESULTS" 50
eme i_on_a1cspline , mat(E2)
eme s_on_a1cspline , mat(E2)
di "`i_on_a1cspline'"
di "`s_on_a1cspline'"