---
author: Rich Jones
generator: pandoc
title: Joint Model
viewport: width=device-width, initial-scale=1.0, user-scalable=yes
---

Rich Jones

July 09 2023

# Overview

This document summarizes analyses conducted for 2022 ΨMCA work group 3,
focusing on the effect of blood sugar (HbA1c) in memory performance
change over 10 years in Black or African-American participants the
Health and Retirement Study (HRS).

The goal of this analysis was to estimate a joint model (a growth curve
model with a survival model) to compare to a model that uses propensity
score weighted analyses (Laura Gibbons is conducting those analyses). To
build up to the joint model, I estimate: a (1) survival model in Stata
to compare to a (2) survival model in Mplus, a (3) growth curve model,
and finally a (4) joint model. Unconditional (i.e., no covariates) and
conditional (i.e., with covariates) models are analyzed for all but the
joint model (only a conditional joint model is estimated). The
covariates used in the adjusted or conditional models are all of those
factors that Laura used in generating the selection weights.

Missing data among the covariates were handled with multiple imputation.
Laura Gibbons performed the imputation and provided the analytic data
set. All analyses in this report report results based on 5 imputations.
The descriptive table uses the non-imputed data.

This analysis uses HRS data, but does not use sampling weights.

## Summary of findings

HbA1c (the spline version) is significantly related to risk of death (HR
1.12 \[1.05, 1.19\] per unit difference in HbA1c spline) in the
unconditional model and this effect is attenuated (HR 1.07) when
controlling for covariates but remains significant, including in the
joint model.

HbA1c spline is not significantly related to *level* (i.e., intercept)
and *slope* (change over 10 years follow-up) in memory performance, and
this is true in the unconditional, conditional, and joint model. The
point estimate suggests lower memory performance with higher HbA1c, and
faster (i.e., more negative) slope with higher HbA1c, but we cannot rule
out that the true effect is no effect given the sample size (confidence
intervals include 0). The effect of HbA1c in intercept is relatively
unchanged between the unconditional, conditional, and joint model. The
effect of HbA1c in slope is about a half of the magnitude in the
unconditional and joint models relative to the unconditional model.

The conditional and joint models are not appreciably different.
Therefore, conditioning the growth model on survival did not influence
how we assess the relationship between HbA1c and memory performance and
change.

Finally, the Cox model estimated with Stata and the continuous time
survival model estimated with Mplus produce slightly different parameter
estimates. I don't have an explanation for that, and therefore all of
the results are somewhat tenuous.

# Contents

1.  [Overview](#overview)
2.  [Analysis](#analysis)
    a.  [Data and processing](#data-and-processing)
    b.  [Survival models](#survival-models)
        i.  [Continuous time survival with multiply imputed data,
            unconditional model in
            Stata](#continuous-time-survival-with-multiply-imputed-data-unconditional-model-in-stata)
        ii. [Continuous time survival with multiply imputed data,
            unconditional model in
            Mplus](#continuous-time-survival-with-multiply-imputed-data-unconditional-model-in-mplus)
    c.  [Latent growth curve model](#latent-growth-curve-model)
        i.  [Unconditional, imputed data](#unconditional-imputed-data)
        ii. [Conditional, imputed data](#conditional-imputed-data)
        iii. [Joint model](#joint-model)
3.  [Summary table](#summary-table)
4.  [Sample descriptive table](#sample-descriptive-table)
5.  [Detailed results and program
    files](#detailed-results-and-program-files)

# Analysis

## Call data

Using LG's file `forsas230518.dta`.

``` {.stata}
. * chose variables for analysis
. use "forsas230518.dta" , clear

. gen _mi_miss=0

. gen _mi_id=0

. mi extract 0, clear

. tempfile foo

. gen m=_Imputation_

. cap drop _mi_miss

. cap drop _mi_id

. cap drop _mi_m 

. * keep if race==2 // Black or African-American (regardless of ethnicity)
. * keep if wg3awave~=. // if observed in our time window
```

``` {.stata}
. * covariates used in generating weights 
. su ///
>    hba1c wg3 gender ///
>    ed3 /// i.ed3 
>    spl_age1-spl_age4 ///
>    cses_index maxparented couple usborn yrenter fullret ///
>    hi25 hinpova rbmi rshlt c095 c103 rcesd  ///
>    smoke3 /// i.smoke3 
>    rhibpe /*rhearte*/ radla riadla rgrossa rmobila rlgmusa ///
>    riadlza rfinea rshltc rdrink rlunge rstroke rcancre ///
>    rarthre rconde routpt rdentst rdrugs rhomcar rhsptim ///
>    rnrstim rhspnit rdoctim rnrsnit rnhmliv hatotb ///
>    hacohort /// i.hacohort 
>    apoe4 ///
>    gender /// i.gender#c.age0 c.age0#c.wg3 i.gender#c.wg3 i.apoe4#c.wg3 i.ed3#c.wg3
>    if missing(hba1c)~=1 & missing(rtr20)~=1

    Variable │        Obs        Mean    Std. dev.       Min        Max
─────────────┼─────────────────────────────────────────────────────────
       hba1c │     36,906    6.166721    1.162206       3.57      14.33
    wg3awave │     36,906    1.949602    1.566848          0          5
      gender │     36,906    1.635181    .4813858          1          2
         ed3 │     36,906    .8405137    .6886004          0          2
    spl_age1 │     36,906    71.06194    8.938769         51        101
─────────────┼─────────────────────────────────────────────────────────
    spl_age2 │     36,906    6.081149    8.093625          0   52.25806
    spl_age3 │     36,906    1.329253    2.531165          0   18.13111
    spl_age4 │     36,906     .282179    .6961501          0   5.494277
  cses_index │     36,906   -.2334339    .8882616  -3.638843   2.164723
 maxparented │     36,183    8.926131    3.554733  -5.397321   18.08955
─────────────┼─────────────────────────────────────────────────────────
      couple │     36,906    .4186311    .4933414          0          1
      usborn │     36,901    .9364245    .2439985          0          1
     yrenter │     36,906    1996.124    4.931397       1992       2008
     fullret │     36,565    .6409408    .4797312          0          1
        hi25 │     36,906    12.77513     53.9055        -25   1544.871
─────────────┼─────────────────────────────────────────────────────────
     hinpova │     36,906    .2378475     .425771          0          1
        rbmi │     36,853    29.61023    6.554338    11.1389       75.8
       rshlt │     36,896    3.176518    1.024689   .6259943   5.979667
        c095 │     36,084    3.186569    .9853053  -.5024406   6.463416
        c103 │     36,092    2.719128    1.047947  -1.310774   6.284186
─────────────┼─────────────────────────────────────────────────────────
       rcesd │     36,902    1.730675    2.074737  -.9217504          8
      smoke3 │     36,873    .7077265    .6949428          0          2
      rhibpe │     36,896    .7885137     .408368          0          1
       radla │     36,899    .4878275    1.062353  -1.975416          5
      riadla │     36,899    .1779733    .5332721  -1.014224          3
─────────────┼─────────────────────────────────────────────────────────
     rgrossa │     36,899     .758243    1.277397  -1.741531          5
     rmobila │     36,898    1.490823    1.615261  -2.341007          5
     rlgmusa │     36,899    1.586925    1.413423  -1.726378          4
     riadlza │     36,899    .3978531    .9516359  -1.462543          5
      rfinea │     36,899    .2856681    .6161086  -1.584405          3
─────────────┼─────────────────────────────────────────────────────────
      rshltc │     36,881    .0477611    .9081867         -4          4
      rdrink │     36,906    .3604292    .4801315          0          1
      rlunge │     36,893    .0898544     .285977          0          1
     rstroke │     36,906    .1216062    .3268349          0          1
     rcancre │     36,881    .1626854    .3690834          0          1
─────────────┼─────────────────────────────────────────────────────────
     rarthre │     36,895    .6863532     .463981          0          1
      rconde │     36,906    2.598277    1.470277          0          8
      routpt │     36,880    .1545282    .3614592          0          1
     rdentst │     36,882    .4285288    .4948722          0          1
      rdrugs │     36,892    .8830912     .321316          0          1
─────────────┼─────────────────────────────────────────────────────────
     rhomcar │     36,874     .128356    .3344903          0          1
     rhsptim │     36,849    .5709714    1.268112  -2.319584         20
     rnrstim │     36,873    .0446255     .284066  -2.581517          7
     rhspnit │     36,813    2.712522     9.99831  -23.43909        270
     rdoctim │     36,315    10.78737    23.93209  -93.29469        900
─────────────┼─────────────────────────────────────────────────────────
     rnrsnit │     36,873     3.78105    42.80817  -131.2104        975
     rnhmliv │      6,151    .0055276    .0741478          0          1
      hatotb │      6,151    159308.1    489753.5    -342160   1.70e+07
    hacohort │      6,151    3.395708    1.084469          0          5
       apoe4 │     36,694    .3713414    .4831701          0          1
─────────────┼─────────────────────────────────────────────────────────
      gender │     36,906    1.635181    .4813858          1          2

. 
. su age0 if wg3awave==0

    Variable │        Obs        Mean    Std. dev.       Min        Max
─────────────┼─────────────────────────────────────────────────────────
        age0 │      8,820    60.84694    8.019928         50         88

. gen cage0=age0-`r(mean)'

. *'
. gen female=gender==2 if missing(gender)~=1

. 
. gen int1 = cage0*female 

. *# None of the wg3awave-interaction covariates are going 
. *# to Rich's analysis because I am only taking baseline covariates
. gen int2 = cage0*wg3awave 

. gen int3 = female*wg3awave 

. gen int4 = apoe4*wg3awave
(212 missing values generated)

. gen int5 = (ed3==0)*wg3awave 

. gen int6 = (ed3==2)*wg3awave 

. 
. gen lths = ed3==0

. gen gths = ed3==2

. 
. local covlist "female lths gths "

. local covlist "`covlist' spl_age1-spl_age4"

. local covlist "`covlist' cses_index maxparented couple usborn"

. center yrenter if wg3awave==0
cyrenter generated as yrenter centered at 1995.846, with 0 missing values, the mean of the group defined by if wg3awave==0

. local covlist "`covlist' cyrenter fullret hi25 hinpova rbmi rshlt c095 c103 rcesd " 

. gen pastsmk = smoke3==1 if missing(smoke3)~=1
(33 missing values generated)

. gen currsmk = smoke3==2 if missing(smoke3)~=1
(33 missing values generated)

. local covlist "`covlist' pastsmk currsmk"

. local covlist "`covlist' rhibpe "

. *** local covlist "`covlist' rhearte" 
. local covlist "`covlist' radla riadla rgrossa rmobila rlgmusa"

. local covlist "`covlist' riadlza rfinea rshltc rdrink rlunge rstroke rcancre" 

. local covlist "`covlist' rarthre rconde routpt rdentst rdrugs rhomcar rhsptim"

. local covlist "`covlist' rnrstim rhspnit rdoctim rnrsnit  hatotb " // rnhmliv <- a constant at baseline 

. 
. /* hacohort is not imputed (always missing when _Imputation_>0)
> gen cohort01 = hacohort==1 | hacohort==0
> gen cohort2 = hacohort==2
> gen cohort4 = hacohort==4
> gen cohort5 = hacohort==5 
> 
> local covlist "`covlist' cohort01 cohort2 cohort4 cohort5"
> */
. 
. 
. local covlist "`covlist' apoe4 int1"

. 
. ** Process data
. 
. cap drop k

. gen k=1

. 
. * forward fill hba1c
. gsort _Imputation_ id wg3awave

. egen iid=group(_Imputation_ id)

. xfill hba1c , i(iid)

. 
. * everdead
. preserve

. keep _Imputation_ id diednext 

. collapse (max) diednext , by(_Imputation_ id)

. rename diednext everdead

. * Everdead is only ever 0 and 1 in _Imputation_ 0
. keep if _Imputation_ == 0
(7,350 observations deleted)

. keep id everdead 

. tempfile everdead 

. save `everdead'
file /var/folders/lq/w3m6z0dj41ngkbbc0204xb7m0000gp/T//S_47104.000003 saved as .dta format

. *'
. restore 

. merge m:1 id using `everdead'

    Result                      Number of obs
    ─────────────────────────────────────────
    Not matched                             0
    Matched                            36,906  (_merge==3)
    ─────────────────────────────────────────

. *'
. table _Imputation_ everdead

──────────┬─────────────
          │    (max)    
_Imputati │   diednext  
on_       │     0      1
──────────┼─────────────
        0 │ 4,889  1,262
        1 │ 4,889  1,262
        2 │ 4,889  1,262
        3 │ 4,889  1,262
        4 │ 4,889  1,262
        5 │ 4,889  1,262
──────────┴─────────────

. 
. 
. tempfile masterdata

. save `masterdata' , replace 
(file /var/folders/lq/w3m6z0dj41ngkbbc0204xb7m0000gp/T//S_47104.000004 not found)
file /var/folders/lq/w3m6z0dj41ngkbbc0204xb7m0000gp/T//S_47104.000004 saved as .dta format

. *'
. * longitudinal data
. keep _Imputation_ id wg3awave rtr20 diednext everdead 

. tempfile longdata_long

. save `longdata_long' , replace 
(file /var/folders/lq/w3m6z0dj41ngkbbc0204xb7m0000gp/T//S_47104.000005 not found)
file /var/folders/lq/w3m6z0dj41ngkbbc0204xb7m0000gp/T//S_47104.000005 saved as .dta format

. *'
. reshape wide  rtr20 diednext , i(_Imputation_ id everdead) j(wg3awave)
(j = 0 1 2 3 4 5)

Data                               Long   ->   Wide
─────────────────────────────────────────────────────────────────────────────
Number of observations           36,906   ->   8,820       
Number of variables                   6   ->   15          
j variable (6 values)          wg3awave   ->   (dropped)
xij variables:
                                  rtr20   ->   rtr200 rtr201 ... rtr205
                               diednext   ->   diednext0 diednext1 ... diednext5
─────────────────────────────────────────────────────────────────────────────

. tempfile longdata_wide 

. save `longdata_wide' 
file /var/folders/lq/w3m6z0dj41ngkbbc0204xb7m0000gp/T//S_47104.000006 saved as .dta format

. *'
. * Covariates 
. use `masterdata' , clear 

. *'
. keep if wg3awave ==0
(28,086 observations deleted)

. keep _Imputation_ id hba1c `covset' 

. *'
. tempfile covariates_baseline 

. save `covariates_baseline'
file /var/folders/lq/w3m6z0dj41ngkbbc0204xb7m0000gp/T//S_47104.000007 saved as .dta format

. * wide analysis file 
. *'
. use `longdata_wide' , clear

. *' 
. merge 1:1 _Imputation_ id using `covariates_baseline' 

    Result                      Number of obs
    ─────────────────────────────────────────
    Not matched                             0
    Matched                             8,820  (_merge==3)
    ─────────────────────────────────────────

. 
. forvalues i=1/6 {
  2. local j=`i'-1
  3. rename diednext`j' dead`i'
  4.    *'
. }

. replace everdead=0 if everdead~=1
(0 real changes made)

. 
. * person years of follow-up
. gen     py = 14 if missing(everdead)~=1

. replace py = 11 if dead6==1
(22 real changes made)

. replace py = 9 if dead5==1
(76 real changes made)

. replace py = 7 if dead4==1
(65 real changes made)

. replace py = 5 if dead3==1
(81 real changes made)

. replace py = 3 if dead2==1
(82 real changes made)

. replace py = 1 if dead1==1
(89 real changes made)

. 
. 
. local knotis = 6.5

. gen a1cspline = 0 if missing(hba1c)~=1

. replace a1cspline = max(0,hba1c-`knotis')
(2,022 real changes made)

. *'
. plot a1cspline hba1c

    7.83 +  
         |                                                                 *
         |                                                                *
         |                                                            *
         |  
         |                                                         *
    a    |                                                     *
    1    |                                                  **
    c    |                                               ***
    s    |                                             **
    p    |                                           ***
    l    |                                        ***
    i    |                                      ***
    n    |                                    **
    e    |                                 ***
         |                              ***
         |                            ***
         |                         ****
         |                       ***
         |                    ****
       0 + *   ****************
          +────────────────────────────────────────────────────────────────+
             3.57               A1C at wg3 baseline                  14.33


. 
. 
. local covset "a1cspline `covlist'"

. 
. * In Mplus, we use a variable TIMECENSORED
. * to indicate whether or not an event occurred.
. * The tricky thing is, time is censored for people
. * who DID NOT experience the event, and it is
. * not censored for people who DID experience the event.
. gen tc=0 // initialize to 0

. replace tc=1 if everdead==0
(6,330 real changes made)

. 
. 
. * mi set 
. 
. tempfile foo

. gen m=_Imputation_

. cap drop _mi_miss

. cap drop _mi_id

. cap drop _mi_m 

. save `foo'
file /var/folders/lq/w3m6z0dj41ngkbbc0204xb7m0000gp/T//S_47104.000008 saved as .dta format

. use `foo'

. mi import flong , m(m) id(id)

. 
. cap drop _merge 

. merge 1:1  id _Imputation_ using covariates.dta 
(label YESNO already defined)

    Result                      Number of obs
    ─────────────────────────────────────────
    Not matched                             0
    Matched                             8,820  (_merge==3)
    ─────────────────────────────────────────

. 
. mi stset py, failure(everdead) 
(regular variables age0 rahispan gender ed3 yrenter hacohort unregistered because not in m=0)
(imputed variables smoke3 rdiabe unregistered because not in m=0)
(448 m=0 obs now marked as incomplete)
(7350 values of regular variable hatotb in m>0 updated to match values in m=0)
(7350 values of regular variable rnhmliv in m>0 updated to match values in m=0)

Survival-time data settings

         Failure event: everdead!=0 & everdead<.
Observed time interval: (0, py]
     Exit on or before: failure

──────────────────────────────────────────────────────────────────────────
      1,470  total observations
          0  exclusions
──────────────────────────────────────────────────────────────────────────
      1,470  observations remaining, representing
        415  failures in single-record/single-failure data
     16,891  total analysis time at risk and under observation
                                                At risk from t =         0
                                     Earliest observed entry t =         0
                                          Last observed exit t =        14

. gen iid=_mi_id

. save usingdata.dta , replace 
file usingdata.dta saved

. 
```

## Survival models

### Continuous time Survival with multiply imputed data, unconditional model in Stata

``` {.stata}
. use usingdata.dta , clear 

. * continuous time survival, multiply imputed data
. mi estimate : stcox a1cspline 

Multiple-imputation estimates                   Imputations       =          5
Cox regression: Breslow method for ties         Number of obs     =      1,470
                                                Average RVI       =     0.0000
                                                Largest FMI       =     0.0000
DF adjustment:   Large sample                   DF:     min       =          .
                                                        avg       =          .
                                                        max       =          .
Model F test:       Equal FMI                   F(   1,      .)   =       8.74
Within VCE type:          OIM                   Prob > F          =     0.0031

─────────────┬────────────────────────────────────────────────────────────────
          _t │ Coefficient  Std. err.      t    P>|t|     [95% conf. interval]
─────────────┼────────────────────────────────────────────────────────────────
   a1cspline │   .1319719   .0446519     2.96   0.003     .0444557    .2194881
─────────────┴────────────────────────────────────────────────────────────────

. *'
. mat B=e(b_mi)

. mat B=B'

. mat list B

symmetric B[1,1]
                  r1
a1cspline  .13197189

. di exp(B[1,1])
1.1410762
```

### Continuous time Survival with multiply imputed data, unconditional model in Mplus

``` {.stata}
. include model1.do 

. * See MPUG ex 6.20
. cd $cwd
/Users/rnj/Library/CloudStorage/Dropbox/FH/FH2022/WG3/PsyMCA-WG3

. use usingdata.dta , clear 

. cap confirm file modelinp1.out

. if _rc~=0 {
.    cap erase inpfiles.dat 
.    rdoc init inpfiles.dat 
.    * run individually but only for purposes of saving inputdatafile 
.    forvalues i=1/5 {
  2.       qui runmplus tc py a1cspline if _Imputation_ == `i' , /// `
>          variable(SURVIVAL = py; TIMECENSORED = tc (0 = NOT 1 = RIGHT);) ///
>            model(py ON a1cspline;) ///
>            log(off) saveinp(model`i') savelog(model`i') ///`
>          saveinputdatafile(inp`i') // `
  3.       r inp`i'.dat //`
  4.    }
.    rdoc close 
.    cap erase modelinp1.inp
.    rdoc init modelinp1.inp
.    r DATA: FILE = inpfiles.dat ; 
.    r TYPE = imputation ;
.    r VARIABLE: NAMES = tc py a1cspline ;
.    r MISSING ARE ALL (-9999) ;
.    r SURVIVAL = py ;
.    r TIMECENSORED = tc (0 = NOT 1 = RIGHT) ;
.    r MODEL: py ON a1cspline ;
.    r 
.    rdoc close 
.    type modelinp1.inp
.    !$mplus_path $cwd/modelinp1.inp
. }

. runmplus , po(modelinp1.out) log(off)
(359 missing values generated)
log suppressed

. mat E1 = r(estimate)

. mat SE1 = r(se)

. mat z1 = r(z)

. runmplus_show_output_segment modelinp1.out "MODEL RESULTS" 9


  MODEL RESULTS
  
                                                      Two-Tailed   Rate of
                      Estimate       S.E.  Est./S.E.    P-Value    Missing
  
   PY         ON
      A1CSPLINE          0.111      0.032      3.500      0.000      0.000
  
  

. eme py_on_a1cspline , mat(E1)
matrix  element py_on_a1cspline column 1 is -> 0.111 (returned as r(r1))

. di exp(`py_on_a1cspline')
1.1173949

. 
```

## Latent growth curve model

### Unconditional, imputed data

``` {.stata}
. include model2.do 

. cap confirm file model2mi.out 

. if _rc~=0 {
.    cap erase inpfiles.dat 
.    rdoc init inpfiles.dat 
.    * run individually but only for purposes of saving inputdatafile 
.    forvalues i=1/5 {
  2.       cap erase inp`i'.dat 
  3.       runmplus rtr200-rtr205 a1cspline if _Imputation_==`i' , est(MLR) ///
>          model(i s | rtr200@0 rtr201@2 rtr202@4 rtr203@6 rtr204@8 rtr205@10; ///
>          i s on a1cspline ;) ///
>          log(off) savelog(model2`i') saveinp(model2`i') saveinputdatafile(inp`i') //`
  4.       r inp`i'.dat //`
  5.    }
.    rdoc close 
.    cap erase model2mi.inp 
.    filefilter model21.inp model2mi.inp , from("FILE = inp1.dat ;") to("FILE = inpfiles.dat ; \n TYPE = imputation ;")
.    type model2mi.inp 
.    !$mplus_path $cwd/model2mi.inp
. }

. runmplus , po(model2mi.out) log(off)
(643 missing values generated)
log suppressed

. mat E2 = r(estimate)

. mat SE2 = r(se)

. mat z2 = r(z)

. runmplus_show_output_segment model2mi.out "MODEL RESULTS" 50


  MODEL RESULTS
  
                                                      Two-Tailed   Rate of
                      Estimate       S.E.  Est./S.E.    P-Value    Missing
  
   I        |
      RTR200             1.000      0.000    999.000    999.000      0.000
      RTR201             1.000      0.000    999.000    999.000      0.000
      RTR202             1.000      0.000    999.000    999.000      0.000
      RTR203             1.000      0.000    999.000    999.000      0.000
      RTR204             1.000      0.000    999.000    999.000      0.000
      RTR205             1.000      0.000    999.000    999.000      0.000
  
   S        |
      RTR200             0.000      0.000    999.000    999.000      0.000
      RTR201             2.000      0.000    999.000    999.000      0.000
      RTR202             4.000      0.000    999.000    999.000      0.000
      RTR203             6.000      0.000    999.000    999.000      0.000
      RTR204             8.000      0.000    999.000    999.000      0.000
      RTR205            10.000      0.000    999.000    999.000      0.000
  
   I        ON
      A1CSPLINE         -0.084      0.079     -1.062      0.288      0.000
  
   S        ON
      A1CSPLINE         -0.017      0.012     -1.333      0.183      0.000
  
   S        WITH
      I                  0.041      0.044      0.927      0.354      0.000
  
   Intercepts
      RTR200             0.000      0.000    999.000    999.000      0.000
      RTR201             0.000      0.000    999.000    999.000      0.000
      RTR202             0.000      0.000    999.000    999.000      0.000
      RTR203             0.000      0.000    999.000    999.000      0.000
      RTR204             0.000      0.000    999.000    999.000      0.000
      RTR205             0.000      0.000    999.000    999.000      0.000
      I                  8.298      0.085     97.505      0.000      0.000
      S                 -0.127      0.011    -11.670      0.000      0.000
  
   Residual Variances
      RTR200             4.243      0.326     13.007      0.000      0.000
      RTR201             7.139      0.385     18.551      0.000      0.000
      RTR202             4.078      0.245     16.658      0.000      0.000
      RTR203             7.208      0.456     15.812      0.000      0.000
      RTR204             4.681      0.330     14.191      0.000      0.000
      RTR205             6.627      0.830      7.981      0.000      0.000
      I                  6.401      0.389     16.473      0.000      0.000
      S                  0.007      0.007      1.064      0.287      0.000
  

. eme i_on_a1cspline , mat(E2)
matrix  element i_on_a1cspline column 1 is -> -0.084 (returned as r(r1))

. eme s_on_a1cspline , mat(E2)
matrix  element s_on_a1cspline column 1 is -> -0.017 (returned as r(r1))

. di "`i_on_a1cspline'"
-0.084

. di "`s_on_a1cspline'"
-0.017

. 
```

### Conditional, imputed data

``` {.stata}
. include model3.do 

. global covset "`covset'"

. * LGC model
. * only have obserations for wg3awave
. cd $cwd
/Users/rnj/Library/CloudStorage/Dropbox/FH/FH2022/WG3/PsyMCA-WG3

. use usingdata.dta , clear 

. keep py tc rtr20* `covset' _Imputation_  //`

. local i=0

. foreach x of varlist $covset {
  2.    if "`x'"~="a1cspline" {
  3.       gen cov`++i' = `x' 
  4.       qui distinct cov`i' //`
  5.       if `r(ndistinct)'>2 { //`
  6.          qui standardize cov`i' , replace //`
  7.       }
  8.       qui distinct cov`i' //`
  9.       if `r(ndistinct)'==2 { //'
 10.          qui center cov`i' , replace //`
 11.       }
 12.       local foo : var lab `x' //`
 13.       drop `x' //`
 14.       la var cov`i' "`x' `foo'"  //
 15.    } 
 16. }
(190 missing values generated)
(1 missing value generated)
(143 missing values generated)
(13 missing values generated)
(1 missing value generated)
(1 missing value generated)
(9 missing values generated)
(9 missing values generated)
(2 missing values generated)
(13 missing values generated)
(2 missing values generated)
(4 missing values generated)
(2 missing values generated)
(3 missing values generated)
(1 missing value generated)
(1 missing value generated)
(2 missing values generated)
(9 missing values generated)
(1 missing value generated)
(14 missing values generated)
(86 missing values generated)
(1 missing value generated)
(60 missing values generated)

. local covrange "cov1-cov`i'" 

. 
. cap confirm file model3mi.out 

. if _rc~=0 {
.    cap erase inpfiles.dat 
.    rdoc init inpfiles.dat 
.    * run individually but only for purposes of saving inputdatafile 
.    forvalues i=1/5 {
  2.       cap erase inp`i'.dat 
  3.       runmplus rtr200-rtr205 a1cspline `covrange' if _Imputation_==`i' , est(MLR) ///
>          model(i s | rtr200@0 rtr201@2 rtr202@4 rtr203@6 rtr204@8 rtr205@10; i s on a1cspline `covrange';) ///
>          log(off) savelog(model3`i') saveinp(model3`i') saveinputdatafile(inp`i') //`
  4.       r inp`i'.dat //`
  5.    }
.    rdoc close 
.    cap erase model3mi.inp 
.    filefilter model31.inp model3mi.inp , from("FILE = inp1.dat ;") to("FILE = inpfiles.dat ; \n TYPE = imputation ;")
.    type model3mi.inp 
.    !$mplus_path $cwd/model3mi.inp
. }

. runmplus , po(model3mi.out) log(off)
(6,491 missing values generated)
log suppressed

. mat E3 = r(estimate)

. mat SE3 = r(se)

. mat z3 = r(z)

. runmplus_show_output_segment model3mi.out "MODEL RESULTS" 150


  MODEL RESULTS
  
                                                      Two-Tailed   Rate of
                      Estimate       S.E.  Est./S.E.    P-Value    Missing
  
   I        |
      RTR200             1.000      0.000    999.000    999.000      0.000
      RTR201             1.000      0.000    999.000    999.000      0.000
      RTR202             1.000      0.000    999.000    999.000      0.000
      RTR203             1.000      0.000    999.000    999.000      0.000
      RTR204             1.000      0.000    999.000    999.000      0.000
      RTR205             1.000      0.000    999.000    999.000      0.000
  
   S        |
      RTR200             0.000      0.000    999.000    999.000      0.000
      RTR201             2.000      0.000    999.000    999.000      0.000
      RTR202             4.000      0.000    999.000    999.000      0.000
      RTR203             6.000      0.000    999.000    999.000      0.000
      RTR204             8.000      0.000    999.000    999.000      0.000
      RTR205            10.000      0.000    999.000    999.000      0.000
  
   I        ON
      A1CSPLINE         -0.091      0.073     -1.256      0.209      0.003
      COV1               0.934      0.148      6.316      0.000      0.004
      COV2              -0.823      0.159     -5.188      0.000      0.016
      COV3               0.657      0.194      3.384      0.001      0.003
      COV4               0.007      0.462      0.015      0.988      0.003
      COV5               0.520      1.542      0.337      0.736      0.001
      COV6              -2.419      2.328     -1.039      0.299      0.001
      COV7               1.725      1.156      1.493      0.136      0.002
      COV8              -0.075      0.094     -0.793      0.428      0.132
      COV9               0.173      0.108      1.599      0.110      0.266
      COV10             -0.021      0.141     -0.148      0.882      0.002
      COV11              0.487      0.249      1.954      0.051      0.008
      COV12              0.102      0.155      0.657      0.511      0.003
      COV13             -0.167      0.167     -1.000      0.318      0.118
      COV14              0.261      0.068      3.851      0.000      0.001
      COV15             -0.618      0.168     -3.685      0.000      0.001
      COV16              0.312      0.073      4.249      0.000      0.010
      COV17             -0.348      0.101     -3.453      0.001      0.009
      COV18              0.082      0.074      1.097      0.273      0.007
      COV19             -0.013      0.069     -0.191      0.848      0.004
      COV20             -0.305      0.081     -3.789      0.000      0.001
      COV21             -0.070      0.144     -0.490      0.624      0.040
      COV22              0.005      0.203      0.023      0.982      0.023
      COV23             -0.105      0.195     -0.539      0.590      0.003
      COV24              0.017      0.191      0.087      0.931      0.003
      COV25             -0.086      0.132     -0.653      0.514      0.005
      COV26             -0.043      0.235     -0.185      0.854      0.003
      COV27             -0.174      0.167     -1.043      0.297      0.003
      COV28              0.169      0.088      1.920      0.055      0.002
      COV29             -0.250      0.151     -1.660      0.097      0.005
      COV30              0.084      0.114      0.735      0.462      0.001
      COV31              0.164      0.078      2.096      0.036      0.015
      COV32              0.207      0.138      1.495      0.135      0.002
      COV33              0.462      0.295      1.567      0.117      0.001
      COV34             -0.435      0.248     -1.749      0.080      0.001
      COV35              0.446      0.233      1.918      0.055      0.003
      COV36             -0.072      0.181     -0.397      0.691      0.001
      COV37              0.069      0.150      0.462      0.644      0.000
      COV38             -0.273      0.172     -1.583      0.113      0.003
      COV39              0.281      0.138      2.042      0.041      0.003
      COV40              0.191      0.208      0.918      0.359      0.004
      COV41              0.047      0.243      0.195      0.845      0.001
      COV42             -0.036      0.079     -0.460      0.646      0.013
      COV43              0.103      0.069      1.496      0.135      0.008
      COV44             -0.085      0.062     -1.384      0.166      0.045
      COV45              0.138      0.058      2.384      0.017      0.097
      COV46             -0.093      0.056     -1.667      0.095      0.019
      COV47              0.093      0.041      2.274      0.023      0.000
      COV48             -0.174      0.136     -1.277      0.202      0.034
      COV49             -0.716      0.347     -2.067      0.039      0.000
      COV50              0.045      0.108      0.417      0.677      0.001
  
   S        ON
      A1CSPLINE         -0.007      0.011     -0.671      0.502      0.001
      COV1              -0.017      0.024     -0.713      0.476      0.001
      COV2              -0.012      0.024     -0.476      0.634      0.014
      COV3               0.015      0.030      0.502      0.616      0.001
      COV4              -0.049      0.072     -0.687      0.492      0.002
      COV5              -0.095      0.231     -0.414      0.679      0.002
      COV6               0.162      0.366      0.444      0.657      0.001
      COV7              -0.080      0.193     -0.415      0.678      0.001
      COV8              -0.011      0.015     -0.716      0.474      0.160
      COV9               0.017      0.017      0.968      0.333      0.277
      COV10              0.019      0.022      0.872      0.383      0.004
      COV11             -0.077      0.041     -1.846      0.065      0.001
      COV12             -0.016      0.024     -0.660      0.509      0.002
      COV13             -0.042      0.027     -1.542      0.123      0.157
      COV14             -0.025      0.009     -2.783      0.005      0.001
      COV15              0.039      0.028      1.404      0.160      0.001
      COV16             -0.007      0.011     -0.629      0.529      0.013
      COV17             -0.006      0.015     -0.395      0.693      0.002
      COV18              0.007      0.012      0.560      0.576      0.001
      COV19             -0.008      0.011     -0.745      0.457      0.003
      COV20              0.005      0.012      0.427      0.669      0.001
      COV21             -0.019      0.022     -0.862      0.389      0.008
      COV22              0.021      0.032      0.664      0.507      0.006
      COV23              0.045      0.031      1.475      0.140      0.001
      COV24              0.047      0.031      1.526      0.127      0.003
      COV25              0.002      0.020      0.080      0.937      0.001
      COV26             -0.043      0.039     -1.102      0.270      0.001
      COV27              0.040      0.028      1.421      0.155      0.002
      COV28             -0.007      0.014     -0.506      0.613      0.001
      COV29             -0.018      0.024     -0.731      0.465      0.002
      COV30             -0.013      0.018     -0.699      0.485      0.002
      COV31              0.011      0.013      0.833      0.405      0.010
      COV32             -0.033      0.021     -1.537      0.124      0.001
      COV33              0.011      0.043      0.265      0.791      0.001
      COV34             -0.007      0.042     -0.160      0.873      0.002
      COV35             -0.036      0.035     -1.034      0.301      0.001
      COV36              0.018      0.029      0.632      0.527      0.002
      COV37             -0.028      0.023     -1.220      0.222      0.003
      COV38              0.026      0.029      0.900      0.368      0.003
      COV39              0.015      0.022      0.656      0.512      0.001
      COV40             -0.014      0.033     -0.429      0.668      0.000
      COV41             -0.045      0.045     -1.001      0.317      0.001
      COV42              0.001      0.013      0.075      0.940      0.003
      COV43              0.014      0.012      1.162      0.245      0.036
      COV44              0.032      0.009      3.747      0.000      0.051
      COV45              0.012      0.010      1.289      0.197      0.073
      COV46             -0.005      0.008     -0.616      0.538      0.072
      COV47              0.010      0.012      0.815      0.415      0.003
      COV48              0.001      0.021      0.026      0.979      0.015
      COV49              0.018      0.053      0.330      0.741      0.002
      COV50             -0.034      0.020     -1.730      0.084      0.001
  
   S        WITH
      I                  0.037      0.034      1.108      0.268      0.001
  
   Intercepts
      RTR200             0.000      0.000    999.000    999.000      0.000
      RTR201             0.000      0.000    999.000    999.000      0.000
      RTR202             0.000      0.000    999.000    999.000      0.000
      RTR203             0.000      0.000    999.000    999.000      0.000
      RTR204             0.000      0.000    999.000    999.000      0.000
      RTR205             0.000      0.000    999.000    999.000      0.000
      I                  7.561      0.410     18.442      0.000      0.007
      S                 -0.086      0.066     -1.315      0.188      0.009
  
   Residual Variances
      RTR200             4.370      0.312     14.025      0.000      0.001
      RTR201             6.933      0.356     19.461      0.000      0.000
      RTR202             4.162      0.237     17.563      0.000      0.000
      RTR203             6.955      0.428     16.239      0.000      0.001
      RTR204             4.753      0.320     14.872      0.000      0.000
      RTR205             6.616      0.752      8.798      0.000      0.001
      I                  2.568      0.263      9.767      0.000      0.001
      S                 -0.004      0.006     -0.624      0.532      0.001
  

. eme i_on_a1cspline , mat(E3)
matrix  element i_on_a1cspline column 1 is -> -0.091 (returned as r(r1))

. eme s_on_a1cspline , mat(E3)
matrix  element s_on_a1cspline column 1 is -> -0.007 (returned as r(r1))

. di "`i_on_a1cspline'"
-0.091

. di "`s_on_a1cspline'"
-0.007

. 
```

### Joint model

``` {.stata}
. include model4.do 

. * Joint model
. cap confirm file model4mi.out

. if _rc~=0 {
.    cap erase inpfiles.dat 
.    rdoc init inpfiles.dat 
.    * run individually but only for purposes of saving inputdatafile 
.    forvalues i=1/5 {
  2.       cap erase inp`i'.dat 
  3.       runmplus tc py rtr200-rtr205 a1cspline `covrange' if _Imputation_==`i' , ///
>         variable(SURVIVAL = py; TIMECENSORED = tc (0 = NOT 1 = RIGHT);) ///
>          model(py on a1cspline `covrange' ; ///
>             i s | rtr200@0 rtr201@2 rtr202@4 rtr203@6 rtr204@8 rtr205@10; ///
>             i s on a1cspline `covrange'; ///
>             py on i s; ) ///
>          log(off) savelog(model4`i') saveinp(model4`i') saveinputdatafile(inp`i') //
  4.       r inp`i'.dat //`
  5.    }
.    rdoc close 
.    cap erase model4mi.inp 
.    filefilter model41.inp model4mi.inp , from("FILE = inp1.dat ;") to("FILE = inpfiles.dat ; \n TYPE = imputation ;")
.    type model4mi.inp 
.    !$mplus_path $cwd/model4mi.inp
. }

. runmplus , po(model4mi.out) log(off)
(6,686 missing values generated)
log suppressed

. mat E4 = r(estimate)

. mat SE4 = r(se)

. mat z4 = r(z)

. runmplus_show_output_segment model4mi.out "MODEL RESULTS" 205


  MODEL RESULTS
  
                                                      Two-Tailed   Rate of
                      Estimate       S.E.  Est./S.E.    P-Value    Missing
  
   I        |
      RTR200             1.000      0.000    999.000    999.000      0.000
      RTR201             1.000      0.000    999.000    999.000      0.000
      RTR202             1.000      0.000    999.000    999.000      0.000
      RTR203             1.000      0.000    999.000    999.000      0.000
      RTR204             1.000      0.000    999.000    999.000      0.000
      RTR205             1.000      0.000    999.000    999.000      0.000
  
   S        |
      RTR200             0.000      0.000    999.000    999.000      0.000
      RTR201             2.000      0.000    999.000    999.000      0.000
      RTR202             4.000      0.000    999.000    999.000      0.000
      RTR203             6.000      0.000    999.000    999.000      0.000
      RTR204             8.000      0.000    999.000    999.000      0.000
      RTR205            10.000      0.000    999.000    999.000      0.000
  
   I          ON
      A1CSPLINE         -0.089      0.073     -1.226      0.220      0.004
      COV1               0.928      0.148      6.271      0.000      0.003
      COV2              -0.823      0.158     -5.197      0.000      0.010
      COV3               0.656      0.194      3.373      0.001      0.004
      COV4              -0.001      0.463     -0.002      0.999      0.002
      COV5               0.524      1.544      0.339      0.734      0.001
      COV6              -2.414      2.329     -1.036      0.300      0.001
      COV7               1.719      1.156      1.487      0.137      0.002
      COV8              -0.075      0.095     -0.796      0.426      0.141
      COV9               0.174      0.108      1.610      0.107      0.268
      COV10             -0.020      0.141     -0.143      0.886      0.002
      COV11              0.481      0.250      1.925      0.054      0.013
      COV12              0.099      0.155      0.643      0.521      0.003
      COV13             -0.161      0.167     -0.964      0.335      0.117
      COV14              0.259      0.068      3.830      0.000      0.001
      COV15             -0.615      0.168     -3.667      0.000      0.001
      COV16              0.313      0.074      4.249      0.000      0.010
      COV17             -0.352      0.101     -3.486      0.000      0.007
      COV18              0.082      0.074      1.096      0.273      0.007
      COV19             -0.013      0.069     -0.185      0.854      0.004
      COV20             -0.304      0.081     -3.769      0.000      0.001
      COV21             -0.068      0.144     -0.472      0.637      0.034
      COV22              0.010      0.202      0.048      0.962      0.019
      COV23             -0.105      0.195     -0.537      0.591      0.003
      COV24              0.017      0.191      0.089      0.929      0.003
      COV25             -0.089      0.132     -0.677      0.498      0.006
      COV26             -0.047      0.235     -0.199      0.842      0.003
      COV27             -0.171      0.167     -1.023      0.306      0.002
      COV28              0.168      0.088      1.905      0.057      0.001
      COV29             -0.247      0.151     -1.640      0.101      0.004
      COV30              0.082      0.114      0.721      0.471      0.001
      COV31              0.167      0.078      2.131      0.033      0.013
      COV32              0.204      0.138      1.473      0.141      0.002
      COV33              0.464      0.295      1.575      0.115      0.001
      COV34             -0.443      0.248     -1.785      0.074      0.001
      COV35              0.441      0.233      1.897      0.058      0.003
      COV36             -0.079      0.181     -0.438      0.662      0.002
      COV37              0.075      0.150      0.500      0.617      0.001
      COV38             -0.271      0.172     -1.572      0.116      0.002
      COV39              0.284      0.138      2.064      0.039      0.002
      COV40              0.187      0.208      0.900      0.368      0.003
      COV41              0.043      0.243      0.177      0.859      0.001
      COV42             -0.038      0.079     -0.477      0.633      0.011
      COV43              0.103      0.068      1.503      0.133      0.006
      COV44             -0.085      0.062     -1.381      0.167      0.035
      COV45              0.138      0.058      2.388      0.017      0.095
      COV46             -0.095      0.055     -1.710      0.087      0.020
      COV47              0.095      0.041      2.338      0.019      0.005
      COV48             -0.180      0.136     -1.321      0.186      0.032
      COV49             -0.713      0.347     -2.055      0.040      0.000
      COV50              0.044      0.108      0.404      0.686      0.001
  
   S          ON
      A1CSPLINE         -0.008      0.011     -0.712      0.476      0.003
      COV1              -0.016      0.024     -0.645      0.519      0.000
      COV2              -0.011      0.024     -0.435      0.664      0.004
      COV3               0.016      0.030      0.524      0.600      0.001
      COV4              -0.047      0.072     -0.657      0.511      0.000
      COV5              -0.095      0.231     -0.411      0.681      0.002
      COV6               0.157      0.367      0.429      0.668      0.002
      COV7              -0.076      0.193     -0.391      0.696      0.001
      COV8              -0.011      0.015     -0.705      0.481      0.175
      COV9               0.017      0.018      0.958      0.338      0.279
      COV10              0.019      0.022      0.867      0.386      0.004
      COV11             -0.074      0.041     -1.794      0.073      0.008
      COV12             -0.015      0.024     -0.634      0.526      0.000
      COV13             -0.043      0.027     -1.589      0.112      0.155
      COV14             -0.025      0.009     -2.759      0.006      0.007
      COV15              0.038      0.028      1.367      0.171      0.000
      COV16             -0.007      0.011     -0.640      0.522      0.004
      COV17             -0.005      0.015     -0.316      0.752      0.000
      COV18              0.007      0.012      0.546      0.585      0.001
      COV19             -0.008      0.011     -0.763      0.445      0.002
      COV20              0.005      0.012      0.399      0.690      0.001
      COV21             -0.020      0.022     -0.894      0.371      0.007
      COV22              0.020      0.032      0.632      0.528      0.003
      COV23              0.046      0.031      1.492      0.136      0.001
      COV24              0.047      0.031      1.531      0.126      0.000
      COV25              0.003      0.021      0.138      0.890      0.005
      COV26             -0.043      0.039     -1.096      0.273      0.001
      COV27              0.039      0.028      1.399      0.162      0.000
      COV28             -0.007      0.014     -0.491      0.624      0.000
      COV29             -0.019      0.025     -0.760      0.447      0.007
      COV30             -0.012      0.018     -0.671      0.502      0.000
      COV31              0.010      0.013      0.773      0.439      0.004
      COV32             -0.032      0.021     -1.509      0.131      0.000
      COV33              0.011      0.043      0.257      0.797      0.002
      COV34             -0.004      0.042     -0.090      0.928      0.000
      COV35             -0.035      0.035     -1.007      0.314      0.004
      COV36              0.020      0.029      0.701      0.483      0.004
      COV37             -0.030      0.023     -1.290      0.197      0.009
      COV38              0.026      0.029      0.882      0.378      0.004
      COV39              0.014      0.022      0.622      0.534      0.000
      COV40             -0.014      0.033     -0.408      0.683      0.000
      COV41             -0.044      0.045     -0.964      0.335      0.000
      COV42              0.001      0.013      0.074      0.941      0.003
      COV43              0.014      0.012      1.196      0.232      0.021
      COV44              0.032      0.009      3.782      0.000      0.022
      COV45              0.012      0.010      1.283      0.200      0.071
      COV46             -0.005      0.008     -0.588      0.557      0.045
      COV47              0.009      0.012      0.759      0.448      0.017
      COV48              0.002      0.021      0.106      0.915      0.012
      COV49              0.016      0.053      0.295      0.768      0.001
      COV50             -0.033      0.020     -1.676      0.094      0.000
  
   PY         ON
      I                 -0.022      0.025     -0.878      0.380      0.002
      S                 -0.278      0.197     -1.416      0.157      0.007
  
   PY         ON
      A1CSPLINE          0.063      0.032      1.957      0.050      0.003
      COV1              -0.335      0.104     -3.214      0.001      0.001
      COV2              -0.012      0.094     -0.123      0.902      0.003
      COV3               0.073      0.135      0.540      0.589      0.002
      COV4              -0.004      0.307     -0.015      0.988      0.006
      COV5              -0.695      1.042     -0.668      0.504      0.000
      COV6               1.768      1.478      1.196      0.232      0.001
      COV7              -1.008      0.691     -1.459      0.144      0.003
      COV8              -0.011      0.056     -0.205      0.838      0.044
      COV9               0.052      0.064      0.808      0.419      0.076
      COV10              0.056      0.096      0.586      0.558      0.001
      COV11              0.313      0.195      1.601      0.109      0.002
      COV12             -0.176      0.094     -1.869      0.062      0.003
      COV13              0.306      0.131      2.332      0.020      0.131
      COV14             -0.304      0.149     -2.046      0.041      0.004
      COV15             -0.103      0.100     -1.028      0.304      0.005
      COV16              0.013      0.048      0.267      0.789      0.004
      COV17              0.225      0.063      3.561      0.000      0.007
      COV18             -0.010      0.042     -0.248      0.804      0.001
      COV19             -0.063      0.041     -1.538      0.124      0.007
      COV20             -0.017      0.047     -0.364      0.716      0.008
      COV21              0.184      0.095      1.941      0.052      0.004
      COV22              0.593      0.116      5.104      0.000      0.006
      COV23             -0.118      0.123     -0.959      0.338      0.008
      COV24             -0.016      0.099     -0.160      0.873      0.001
      COV25             -0.086      0.057     -1.489      0.136      0.001
      COV26             -0.121      0.126     -0.961      0.337      0.002
      COV27              0.217      0.100      2.182      0.029      0.002
      COV28             -0.203      0.057     -3.553      0.000      0.001
      COV29              0.107      0.070      1.521      0.128      0.002
      COV30              0.025      0.055      0.452      0.651      0.002
      COV31             -0.130      0.049     -2.673      0.008      0.017
      COV32             -0.050      0.090     -0.549      0.583      0.001
      COV33             -0.026      0.137     -0.186      0.853      0.005
      COV34             -0.066      0.121     -0.550      0.583      0.013
      COV35              0.080      0.120      0.666      0.505      0.004
      COV36             -0.267      0.114     -2.341      0.019      0.002
      COV37              0.272      0.077      3.529      0.000      0.013
      COV38             -0.116      0.120     -0.970      0.332      0.005
      COV39             -0.169      0.091     -1.849      0.064      0.002
      COV40             -0.119      0.161     -0.742      0.458      0.004
      COV41              0.128      0.106      1.212      0.226      0.021
      COV42              0.065      0.037      1.740      0.082      0.018
      COV43              0.051      0.029      1.786      0.074      0.004
      COV44              0.056      0.027      2.072      0.038      0.016
      COV45              0.043      0.025      1.700      0.089      0.143
      COV46             -0.097      0.056     -1.737      0.082      0.000
      COV47             -0.025      0.050     -0.504      0.614      0.003
      COV48              0.006      0.086      0.073      0.942      0.125
      COV49             -0.005      0.201     -0.024      0.981      0.004
      COV50              0.080      0.061      1.320      0.187      0.001
  
   Intercepts
      RTR200             0.000      0.000    999.000    999.000      0.000
      RTR201             0.000      0.000    999.000    999.000      0.000
      RTR202             0.000      0.000    999.000    999.000      0.000
      RTR203             0.000      0.000    999.000    999.000      0.000
      RTR204             0.000      0.000    999.000    999.000      0.000
      RTR205             0.000      0.000    999.000    999.000      0.000
      I                  7.573      0.410     18.464      0.000      0.008
      S                 -0.091      0.066     -1.389      0.165      0.008
  
   Residual Variances
      RTR200             4.190      0.230     18.205      0.000      0.001
      RTR201             6.915      0.356     19.404      0.000      0.000
      RTR202             4.170      0.237     17.572      0.000      0.000
      RTR203             6.958      0.429     16.236      0.000      0.001
      RTR204             4.695      0.316     14.868      0.000      0.001
      RTR205             6.457      0.717      9.002      0.000      0.001
      I                  2.769      0.167     16.542      0.000      0.001
      S                  0.002      0.004      0.517      0.605      0.004
  
  

. 
```

# Summary Table

  --------------------------------------------------------------------------------------------------------------------------------
  Model                                Unconditional                 Conditional                 Joint               
                                                                                                 model,conditional   
  ------------------------------------ --------------- ------------- ------------- ------------- ------------------- -------------
  **Parameter**                        **Est (95%      **p-value**   **Est (95%    **p-value**   **Est (95% CI)**    **p-value**
                                       CI)**                         CI)**                                           

  ***Growth curve parameters***                                                                                      

  Intercept                            8.298 (8.131,   \<.001        7.561 (6.757, \<.001        7.573 (6.769,       \<.001
                                       8.465)                        8.365)                      8.377)              

  Slope (conditional mean)             -0.127 (-0.149, \<.001        -0.086        .19           -0.091 (-0.220,     .16
                                       -0.105)                       (-0.215,                    0.038)              
                                                                     0.043)                                          

  ***HbA1c effects in growth curve***                                                                                

  HbA1c spline in intercept            -0.084 (-0.239, .29           -0.091        .21           -0.089 (-0.232,     .22
                                       0.071)                        (-0.234,                    0.054)              
                                                                     0.052)                                          

  HbA1c spline in slope                -0.017 (-0.041, .18           -0.007        .50           -0.008 (-0.030,     .48
                                       0.007)                        (-0.029,                    0.014)              
                                                                     0.015)                                          

  ***HbA1c effects in survival***                                                                                    

  Hazard ratio per unit increase in    1.117 (1.049,   \<.001        1.069 (1.002, .04           1.065 (1.000,       .050
  HbA1c spline                         1.190)                        1.141)                      1.134)              
  --------------------------------------------------------------------------------------------------------------------------------

# Summary Table

Tables of participant characteristics are available as PDFs, which can
be read into modern versions of Word and edited.

``` {.stata}
  2. table1 `set`j'are' , iqr(`set`j'iqr') missing inashell title("Table `j'. `set`j'is'") file(table`j')
  3. }

Table 1. Participant sociodemographics

                                                                                       Mean  (SD)       Observed  
Characteristic                                                                         or n  or (%)     range     
────────────────────────────────────────────────────────────────────────────────────────────────────────────────
                                                                                                                  
Total [n (%)]                                                                          1470  (100)                
                                                                                                                  
Spl_Age1 [M (SD)]                                                                      68.4  (9.4)      [51.0-96.0]
                                                                                                                  
Female [n (%)]                                                                                                    
0                                                                                       552  (37.6)               
1                                                                                       918  (62.4)               
                                                                                                                  
Lths [n (%)]                                                                                                      
0                                                                                       941  (64.0)               
1                                                                                       529  (36.0)               
                                                                                                                  
Gths [n (%)]                                                                                                      
0                                                                                      1234  (83.9)               
1                                                                                       236  (16.1)               
                                                                                                                  
cses measures:Childhood SES index [M (SD)]                                             -0.2  (0.9)      [-3.6-2.2]
                                                                                                                  
inpova: whether in poverty-w inst       [n (%)]                                                                   
0 0.hh inc above pov thresh                                                            1140  (77.6)               
1 1.hh inc below pov thresh                                                             330  (22.4)               
                                                                                                                  
atotb: total all assets inc. 2nd hm--cross-wave         [Median (P'tile 25-75)]     51600.0  [1800.0-159800.0] [-1.4e+05-1.5e+0

hitot minus 25 000  divided by 1000 [Median (P'tile 25-75)]                            -3.0  [-13.6-20.6] [-25.0-1544.9]
                                                                                                                  
Completely retired  from rsayret [n (%)]                                                                          
0 Not Part                                                                              578  (43.6)               
1 Completely                                                                            749  (56.4)               
                                                                                                                  
Maximum of rameduc and rafeduc [M (SD)]                                                 9.0  (3.6)      [0.0-17.0]
                                                                                                                  
Whether coupled or partnered [n (%)]                                                                              
0 0.no                                                                                  795  (54.1)               
1 1.yes                                                                                 675  (45.9)               
                                                                                                                  
BORN IN THE U.S. [n (%)]                                                                                          
0 0.no                                                                                   95  (6.5)                
1 1.yes                                                                                1374  (93.5)               
                                                                                                                  
yrenter centered at 1995.846 [M (SD)]                                                   0.0  (4.8)      [-3.8-12.2]
                                                                                                                  
────────────────────────────────────────────────────────────────────────────────────────────────────────────────
Notes:  143 missing on fullret, 190 missing on maxparented, 1 missing on usborn. 

table1 results are stored in mata workspace in matrix 'results' , and in 

         table1.txt as tab delimited
         table1.csv as CSV
         table1.tex as LaTeX
         table1.pdf as PDF from LaTeX

Table 2. Participant personal characteristics

                                                    Mean  (SD)       Observed  
Characteristic                                      or n  or (%)     range     
─────────────────────────────────────────────────────────────────────────────
                                                                               
Total [n (%)]                                       1470  (100)                
                                                                               
self-reported body mass index=kg m2 [M (SD)]        29.7  (6.6)      [12.3-70.9]
                                                                               
Any APOE-E4 alleles [n (%)]                                                    
0 0.no                                               879  (62.3)               
1 1.yes                                              531  (37.7)               
                                                                               
─────────────────────────────────────────────────────────────────────────────
Notes:  13 missing on rbmi, 60 missing on apoe4. 

table1 results are stored in mata workspace in matrix 'results' , and in 

         table2.txt as tab delimited
         table2.csv as CSV
         table2.tex as LaTeX
         table2.pdf as PDF from LaTeX

Table 3. Participant self-reported health and functioning

                                                               Mean  (SD)       Observed  
Characteristic                                                 or n  or (%)     range     
────────────────────────────────────────────────────────────────────────────────────────
                                                                                          
Total [n (%)]                                                  1470  (100)                
                                                                                          
self-report of health [n (%)]                                                             
1 1.excellent                                                    71  (4.8)                
2 2.very good                                                   332  (22.6)               
3 3.good                                                        469  (31.9)               
4 4.fair                                                        432  (29.4)               
5 5.poor                                                        165  (11.2)               
                                                                                          
cesd score [n (%)]                                                                        
0                                                               518  (35.3)               
1                                                               381  (25.9)               
2                                                               185  (12.6)               
3                                                               100  (6.8)                
4                                                                69  (4.7)                
5                                                                75  (5.1)                
6                                                                69  (4.7)                
7                                                                46  (3.1)                
8                                                                26  (1.8)                
                                                                                          
some diff-adls  0-5 [n (%)]                                                               
0                                                              1148  (78.1)               
1                                                               161  (11.0)               
2                                                                68  (4.6)                
3                                                                46  (3.1)                
4                                                                35  (2.4)                
5                                                                12  (0.8)                
                                                                                          
some diff-iadls: w2 onwards  0-3 [n (%)]                                                  
0                                                              1334  (90.7)               
1                                                                97  (6.6)                
2                                                                20  (1.4)                
3                                                                19  (1.3)                
                                                                                          
iadlza: some diff-iadls: w2 onwards  0-5        [n (%)]                                   
0                                                              1218  (82.9)               
1                                                               129  (8.8)                
2                                                                64  (4.4)                
3                                                                32  (2.2)                
4                                                                12  (0.8)                
5                                                                15  (1.0)                
                                                                                          
grossa: walk1 r clim1 bed bath 0-5      [n (%)]                                           
0                                                              1012  (68.8)               
1                                                               188  (12.8)               
2                                                               112  (7.6)                
3                                                                80  (5.4)                
4                                                                49  (3.3)                
5                                                                29  (2.0)                
                                                                                          
mobila: some diff-mobility  0-5         [n (%)]                                           
0                                                               645  (43.9)               
1                                                               296  (20.1)               
2                                                               196  (13.3)               
3                                                               126  (8.6)                
4                                                               120  (8.2)                
5                                                                87  (5.9)                
                                                                                          
lgmusa: some diff-large muscle  0-4     [n (%)]                                           
0                                                               521  (35.4)               
1                                                               294  (20.0)               
2                                                               243  (16.5)               
3                                                               246  (16.7)               
4                                                               166  (11.3)               
                                                                                          
finea: dime eat dress  0-3      [n (%)]                                                   
0                                                              1202  (81.8)               
1                                                               192  (13.1)               
2                                                                65  (4.4)                
3                                                                11  (0.7)                
                                                                                          
shltc: change in self-reported hlth [n (%)]                                               
-4                                                                1  (0.1)                
-3                                                               11  (0.8)                
-2                                                               39  (2.7)                
-1                                                              260  (17.8)               
0                                                               769  (52.8)               
1                                                               303  (20.8)               
2                                                                62  (4.3)                
3                                                                11  (0.8)                
4                                                                 1  (0.1)                
                                                                                          
────────────────────────────────────────────────────────────────────────────────────────
Notes:  1 missing on rshlt, 1 missing on rcesd, 13 missing on rshltc. 

table1 results are stored in mata workspace in matrix 'results' , and in 

         table3.txt as tab delimited
         table3.csv as CSV
         table3.tex as LaTeX
         table3.pdf as PDF from LaTeX

Table 4. Participant health and disease-related characteristics

                                                 Mean  (SD)       Observed  
Characteristic                                   or n  or (%)     range     
──────────────────────────────────────────────────────────────────────────
                                                                            
Total [n (%)]                                    1470  (100)                
                                                                            
r ever had high blood pressure [n (%)]                                      
0 0.no                                            388  (26.4)               
1 1.yes                                          1080  (73.6)               
                                                                            
lunge: r ever had lung disease [n (%)]                                      
0 0.no                                           1366  (93.1)               
1 1.yes                                           102  (6.9)                
                                                                            
stroke: r ever had stroke [n (%)]                                           
0 0.no                                           1315  (89.5)               
1 1.yes                                           155  (10.5)               
                                                                            
Rcancre [n (%)]                                                             
0 0.no                                           1282  (87.4)               
1 1.yes                                           184  (12.6)               
                                                                            
arthre: r ever had arthritis [n (%)]                                        
0 0.no                                            551  (37.5)               
1 1.yes                                           917  (62.5)               
                                                                            
conde: sum of conditions ever had [n (%)]                                   
0                                                 142  (9.7)                
1                                                 300  (20.4)               
2                                                 402  (27.3)               
3                                                 345  (23.5)               
4                                                 163  (11.1)               
5                                                  81  (5.5)                
6                                                  30  (2.0)                
7                                                   7  (0.5)                
                                                                            
──────────────────────────────────────────────────────────────────────────
Notes:  2 missing on rhibpe, 2 missing on rlunge, 4 missing on rcancre, 2 missing on rarthre. 

table1 results are stored in mata workspace in matrix 'results' , and in 

         table4.txt as tab delimited
         table4.csv as CSV
         table4.tex as LaTeX
         table4.pdf as PDF from LaTeX

Table 5. Participant treatments and services utilization characteristics

                                                                                 Mean  (SD)       Observed  
Characteristic                                                                   or n  or (%)     range     
──────────────────────────────────────────────────────────────────────────────────────────────────────────
                                                                                                            
Total [n (%)]                                                                    1470  (100)                
                                                                                                            
outpt: outpatient surgry  prv 2 yrs     [n (%)]                                                             
0 0.no                                                                           1233  (84.0)               
1 1.yes                                                                           234  (16.0)               
                                                                                                            
dentst: dental visit  prv 2 yrs         [n (%)]                                                             
0 0.no                                                                            852  (58.0)               
1 1.yes                                                                           617  (42.0)               
                                                                                                            
drugs: reg take rx  prv 2 yrs   [n (%)]                                                                     
0 0.no                                                                            218  (14.8)               
1 1.yes                                                                          1251  (85.2)               
                                                                                                            
homcar: home hlth care  prv 2 yrs       [n (%)]                                                             
0 0.no                                                                           1317  (89.7)               
1 1.yes                                                                           151  (10.3)               
                                                                                                            
hsptim: No. hospital stays  prv 2 yrs     [Median (P'tile 25-75)]                 0.0  [0.0-1.0]  [0.0-7.0] 
                                                                                                            
nrstim: No. nurs home stays  prv 2 yrs    [n (%)]                                                           
0                                                                                1449  (98.6)               
1                                                                                  16  (1.1)                
2                                                                                   3  (0.2)                
3                                                                                   1  (0.1)                
                                                                                                            
hspnit: No. nights in hosp  prv 2 yrs     [Median (P'tile 25-75)]                 0.0  [0.0-1.0]  [0.0-180.0]
                                                                                                            
doctim: No. doctor vists  prv 2 yrs       [Median (P'tile 25-75)]                 7.0  [3.0-12.0] [0.0-350.0]
                                                                                                            
nrsnit: No. nights in nurs home  prv 2 yrs        [Median (P'tile 25-75)]         0.0  [0.0-0.0]  [0.0-180.0]
                                                                                                            
──────────────────────────────────────────────────────────────────────────────────────────────────────────
Notes:  3 missing on routpt, 1 missing on rdentst, 1 missing on rdrugs, 2 missing on rhomcar, 9 missing on rhsptim, 1 missing o

table1 results are stored in mata workspace in matrix 'results' , and in 

         table5.txt as tab delimited
         table5.csv as CSV
         table5.tex as LaTeX
         table5.pdf as PDF from LaTeX

Table 6. Participant reported sensory impairment

                             Mean  (SD)       Observed  
Characteristic               or n  or (%)     range     
──────────────────────────────────────────────────────
                                                        
Total [n (%)]                1470  (100)                
                                                        
RATE EYESIGHT [n (%)]                                   
1 1                            84  (5.7)                
2 2                           261  (17.8)               
3 3                           642  (43.7)               
4 4                           352  (23.9)               
5 5                           124  (8.4)                
6 6                             7  (0.5)                
                                                        
RATE HEARING [n (%)]                                    
1 1                           252  (17.1)               
2 2                           377  (25.6)               
3 3                           556  (37.8)               
4 4                           224  (15.2)               
5 5                            61  (4.1)                
                                                        
──────────────────────────────────────────────────────

table1 results are stored in mata workspace in matrix 'results' , and in 

         table6.txt as tab delimited
         table6.csv as CSV
         table6.tex as LaTeX
         table6.pdf as PDF from LaTeX

Table 7. Participant self-reported health behaviors

                                                Mean  (SD)       Observed  
Characteristic                                  or n  or (%)     range     
─────────────────────────────────────────────────────────────────────────
                                                                           
Total [n (%)]                                   1470  (100)                
                                                                           
Pastsmk [n (%)]                                                            
0                                                869  (59.5)               
1                                                592  (40.5)               
                                                                           
Currsmk [n (%)]                                                            
0                                               1222  (83.6)               
1                                                239  (16.4)               
                                                                           
drink: r ever drinks any alcohol [n (%)]                                   
0 0.no                                           950  (64.6)               
1 1.yes                                          520  (35.4)               
                                                                           
─────────────────────────────────────────────────────────────────────────
Notes:  9 missing on pastsmk, 9 missing on currsmk. 

table1 results are stored in mata workspace in matrix 'results' , and in 

         table7.txt as tab delimited
         table7.csv as CSV
         table7.tex as LaTeX
         table7.pdf as PDF from LaTeX
```

# Detailed results and program files

## PDF Tables on GitHub

[Table 1 Participant
sociodemographics](https://github.com/rnj0nes/PsiMCA22WG3/blob/main/table1.pdf)\
[Table 2 Participant personal
characteristics](https://github.com/rnj0nes/PsiMCA22WG3/blob/main/table2.pdf)\
[Table 3 Participant self-reported health and
functioning](https://github.com/rnj0nes/PsiMCA22WG3/blob/main/table3.pdf)\
[Table 4 Participant health and disease-related
characteristics](https://github.com/rnj0nes/PsiMCA22WG3/blob/main/table4.pdf)\
[Table 5 Participant treatments and services utilization
characteristics](https://github.com/rnj0nes/PsiMCA22WG3/blob/main/table5.pdf)\
[Table 6 Participant reported sensory
impairment](https://github.com/rnj0nes/PsiMCA22WG3/blob/main/table6.pdf)\
[Table 7 Participant self-reported health
behaviors](https://github.com/rnj0nes/PsiMCA22WG3/blob/main/table7.pdf)

## Code files on GitHub

### Master

[This domd
file](https://github.com/rnj0nes/PsiMCA22WG3/blob/main/jm-03-continuous_time_survival-LG-data-20230706.domd)\
[This domd file, HTML
output](https://github.com/rnj0nes/PsiMCA22WG3/blob/main/jm-03-continuous_time_survival-LG-data-20230706.html)

### Continuous time survival, unconditional

[Model 1 inp
file](https://github.com/rnj0nes/PsiMCA22WG3/blob/main/modelinp1.inp)\
[Model 1 out
file](https://github.com/rnj0nes/PsiMCA22WG3/blob/main/modelinp1.out)

### Continuous time survival, conditional

[Model 5 do file (runmplus
codes)](https://github.com/rnj0nes/PsiMCA22WG3/blob/main/model5.do)\
[Model 5 inp
file](https://github.com/rnj0nes/PsiMCA22WG3/blob/main/model5mi.inp)\
[Model 5 out
file](https://github.com/rnj0nes/PsiMCA22WG3/blob/main/model5mi.out)

### Latent growth curve model, unconditional

[Model 2 do file (runmplus
codes)](https://github.com/rnj0nes/PsiMCA22WG3/blob/main/model2.do)\
[Model 2 inp
file](https://github.com/rnj0nes/PsiMCA22WG3/blob/main/model2mi.inp)\
[Model 2 out
file](https://github.com/rnj0nes/PsiMCA22WG3/blob/main/model2mi.out)

### Latent growth curve model, conditional

[Model 3 do file (runmplus
codes)](https://github.com/rnj0nes/PsiMCA22WG3/blob/main/model3mi.do)\
[Model 3 inp
file](https://github.com/rnj0nes/PsiMCA22WG3/blob/main/model3mi.inp)\
[Model 3 out
file](https://github.com/rnj0nes/PsiMCA22WG3/blob/main/model3mi.out)

### Joint model

[Model 4 do file (runmplus
codes)](https://github.com/rnj0nes/PsiMCA22WG3/blob/main/m%3Cbr%3Eodel4mi.do)
[Model 4 inp
file](https://github.com/rnj0nes/PsiMCA22WG3/blob/main/model4mi.inp)\
[Model 4 out
file](https://github.com/rnj0nes/PsiMCA22WG3/blob/main/model4mi.out)

(fin)
