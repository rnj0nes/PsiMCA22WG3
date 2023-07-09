* See MPUG ex 6.20
cd $cwd
use usingdata.dta , clear 
cap confirm file modelinp1.out
if _rc~=0 {
   cap erase inpfiles.dat 
   rdoc init inpfiles.dat 
   * run individually but only for purposes of saving inputdatafile 
   forvalues i=1/5 {
      qui runmplus tc py a1cspline if _Imputation_ == `i' , /// `
         variable(SURVIVAL = py; TIMECENSORED = tc (0 = NOT 1 = RIGHT);) ///
   	   model(py ON a1cspline;) ///
   	   log(off) saveinp(model`i') savelog(model`i') ///`
         saveinputdatafile(inp`i') // `
      r inp`i'.dat //`
   }
   rdoc close 
   cap erase modelinp1.inp
   rdoc init modelinp1.inp
   r DATA: FILE = inpfiles.dat ; 
   r TYPE = imputation ;
   r VARIABLE: NAMES = tc py a1cspline ;
   r MISSING ARE ALL (-9999) ;
   r SURVIVAL = py ;
   r TIMECENSORED = tc (0 = NOT 1 = RIGHT) ;
   r MODEL: py ON a1cspline ;
   r 
   rdoc close 
   type modelinp1.inp
   !$mplus_path $cwd/modelinp1.inp
}
runmplus , po(modelinp1.out) log(off)
mat E1 = r(estimate)
mat SE1 = r(se)
mat z1 = r(z)
runmplus_show_output_segment modelinp1.out "MODEL RESULTS" 9
eme py_on_a1cspline , mat(E1)
di exp(`py_on_a1cspline')
