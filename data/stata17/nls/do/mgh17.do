/* NIST/ITL StRD
Dataset Name:  MGH17             (MGH17.dat)

File Format:   ASCII
               Starting Values   (lines 41 to 45)
               Certified Values  (lines 41 to 50)
               Data              (lines 61 to 93)

Procedure:     Nonlinear Least Squares Regression

Description:   This problem was found to be difficult for some very
               good algorithms.

               See More, J. J., Garbow, B. S., and Hillstrom, K. E.
               (1981).  Testing unconstrained optimization software.
               ACM Transactions on Mathematical Software. 7(1):
               pp. 17-41.

Reference:     Osborne, M. R. (1972).  
               Some aspects of nonlinear least squares 
               calculations.  In Numerical Methods for Nonlinear 
               Optimization, Lootsma (Ed).  
               New York, NY:  Academic Press, pp. 171-189.
 
Data:          1 Response  (y)
               1 Predictor (x)
               33 Observations
               Average Level of Difficulty
               Generated Data

Model:         Exponential Class
               5 Parameters (b1 to b5)

               y = b1 + b2*exp[-x*b4] + b3*exp[-x*b5]  +  e



          Starting values                  Certified Values

        Start 1     Start 2           Parameter     Standard Deviation
  b1 =     50         0.5          3.7541005211E-01  2.0723153551E-03
  b2 =    150         1.5          1.9358469127E+00  2.2031669222E-01
  b3 =   -100        -1           -1.4646871366E+00  2.2175707739E-01
  b4 =      1          0.01        1.2867534640E-02  4.4861358114E-04
  b5 =      2          0.02        2.2122699662E-02  8.9471996575E-04

Residual Sum of Squares:                    5.4648946975E-05
Residual Standard Deviation:                1.3970497866E-03
Degrees of Freedom:                                28
Number of Observations:                            33
*/

clear
program drop _all

scalar N         = 33
scalar df_r      = 28
scalar df_m      = 5

scalar rss       = 5.4648946975E-05
scalar rmse      = 1.3970497866E-03

scalar b1        = 3.7541005211E-01  
scalar seb1      = 2.0723153551E-03
scalar b2        = 1.9358469127E+00  
scalar seb2      = 2.2031669222E-01
scalar b3        = -1.4646871366E+00  
scalar seb3      = 2.2175707739E-01
scalar b4        = 1.2867534640E-02  
scalar seb4      = 4.4861358114E-04
scalar b5        = 2.2122699662E-02  
scalar seb5      = 8.9471996575E-04

program define nlbprog
        version 6
        if "`1'"=="?" {
                global S_1 "b1 b2 b3 b4 b5"

* first set of initial values do not lead to convergence for original nl
                global b1 .5
                global b2 1.5         
                global b3 -1         
                global b4 .01      
                global b5 .02     
                exit
        }

        replace `1' = $b1 + $b2*exp(-(x*$b4)) + $b3*exp(-(x*$b5))
end

qui input double(y x)
      8.440000E-01    0.000000E+00
      9.080000E-01    1.000000E+01
      9.320000E-01    2.000000E+01
      9.360000E-01    3.000000E+01
      9.250000E-01    4.000000E+01
      9.080000E-01    5.000000E+01
      8.810000E-01    6.000000E+01
      8.500000E-01    7.000000E+01
      8.180000E-01    8.000000E+01
      7.840000E-01    9.000000E+01
      7.510000E-01    1.000000E+02
      7.180000E-01    1.100000E+02
      6.850000E-01    1.200000E+02
      6.580000E-01    1.300000E+02
      6.280000E-01    1.400000E+02
      6.030000E-01    1.500000E+02
      5.800000E-01    1.600000E+02
      5.580000E-01    1.700000E+02
      5.380000E-01    1.800000E+02
      5.220000E-01    1.900000E+02
      5.060000E-01    2.000000E+02
      4.900000E-01    2.100000E+02
      4.780000E-01    2.200000E+02
      4.670000E-01    2.300000E+02
      4.570000E-01    2.400000E+02
      4.480000E-01    2.500000E+02
      4.380000E-01    2.600000E+02
      4.310000E-01    2.700000E+02
      4.240000E-01    2.800000E+02
      4.200000E-01    2.900000E+02
      4.140000E-01    3.000000E+02
      4.110000E-01    3.100000E+02
      4.060000E-01    3.200000E+02
end

nl ( y = {b1} + {b2}*exp(-x*{b4}) + {b3}*exp(-x*{b5}) ), ///	
	init(b1 .5 b2 1.5 b3 -1 b4 0.01 b5 0.02) eps(1e-10) 

assert N    == e(N)
assert df_r == e(df_r)
assert df_m == e(df_m)

lrecomp [b1]_b[_cons] b1 [b2]_b[_cons] b2 [b3]_b[_cons] b3 /*
*/ [b4]_b[_cons] b4 [b5]_b[_cons] b5 () /*
*/ [b1]_se[_cons] seb1 [b2]_se[_cons] seb2 [b3]_se[_cons] seb3 /*
*/ [b4]_se[_cons] seb4 [b5]_se[_cons] seb5 () /*
*/ e(rmse) rmse e(rss) rss


