global covset "`covset'"
* LGC model
* only have obserations for wg3awave
cd $cwd
use usingdata.dta , clear 
keep py tc rtr20* `covset' _Imputation_  //`
local i=0
foreach x of varlist $covset {
   if "`x'"~="a1cspline" {
      gen cov`++i' = `x' 
      qui distinct cov`i' //`
      if `r(ndistinct)'>2 { //`
         qui standardize cov`i' , replace //`
      }
      qui distinct cov`i' //`
      if `r(ndistinct)'==2 { //'
         qui center cov`i' , replace //`
      }
      local foo : var lab `x' //`
      drop `x' //`
      la var cov`i' "`x' `foo'"  //
   } 
}
local covrange "cov1-cov`i'" 

cap confirm file model3mi.out 
if _rc~=0 {
   cap erase inpfiles.dat 
   rdoc init inpfiles.dat 
   * run individually but only for purposes of saving inputdatafile 
   forvalues i=1/5 {
      cap erase inp`i'.dat 
      runmplus rtr200-rtr205 a1cspline `covrange' if _Imputation_==`i' , est(MLR) ///
         model(i s | rtr200@0 rtr201@2 rtr202@4 rtr203@6 rtr204@8 rtr205@10; i s on a1cspline `covrange';) ///
         log(off) savelog(model3`i') saveinp(model3`i') saveinputdatafile(inp`i') //`
      r inp`i'.dat //`
   }
   rdoc close 
   cap erase model3mi.inp 
   filefilter model31.inp model3mi.inp , from("FILE = inp1.dat ;") to("FILE = inpfiles.dat ; \n TYPE = imputation ;")
   type model3mi.inp 
   !$mplus_path $cwd/model3mi.inp
}
runmplus , po(model3mi.out) log(off)
mat E3 = r(estimate)
mat SE3 = r(se)
mat z3 = r(z)
runmplus_show_output_segment model3mi.out "MODEL RESULTS" 150
eme i_on_a1cspline , mat(E3)
eme s_on_a1cspline , mat(E3)
di "`i_on_a1cspline'"
di "`s_on_a1cspline'"