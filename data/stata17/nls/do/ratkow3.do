/* NIST/ITL StRD
Dataset Name:  Ratkowsky3        (Ratkowsky3.dat)

File Format:   ASCII
               Starting Values   (lines 41 to 44)
               Certified Values  (lines 41 to 49)
               Data              (lines 61 to 75)

Procedure:     Nonlinear Least Squares Regression

Description:   This model and data are an example of fitting  
               sigmoidal growth curves taken from Ratkowsky (1983).  
               The response variable is the dry weight of onion bulbs 
               and tops, and the predictor variable is growing time. 


Reference:     Ratkowsky, D.A. (1983).  
               Nonlinear Regression Modeling.
               New York, NY:  Marcel Dekker, pp. 62 and 88.





Data:          1 Response  (y = onion bulb dry weight)
               1 Predictor (x = growing time)
               15 Observations
               Higher Level of Difficulty
               Observed Data

Model:         Exponential Class
               4 Parameters (b1 to b4)

               y = b1 / ((1+exp[b2-b3*x])**(1/b4))  +  e



          Starting Values                  Certified Values
 
        Start 1     Start 2           Parameter     Standard Deviation
  b1 =   100         700           6.9964151270E+02  1.6302297817E+01
  b2 =    10           5           5.2771253025E+00  2.0828735829E+00
  b3 =     1           0.75        7.5962938329E-01  1.9566123451E-01
  b4 =     1           1.3         1.2792483859E+00  6.8761936385E-01
 
Residual Sum of Squares:                    8.7864049080E+03
Residual Standard Deviation:                2.8262414662E+01
Degrees of Freedom:                                9  *** should be 11 ***
Number of Observations:                           15 
*/ 
 
clear
program drop _all

scalar N         = 15 
scalar df_r      = 11
scalar df_m      = 4

scalar rss       = 8.7864049080E+03
scalar rmse      = 2.8262414662E+01

scalar b1        = 6.9964151270E+02  
scalar seb1      = 1.6302297817E+01
scalar b2        = 5.2771253025E+00  
scalar seb2      = 2.0828735829E+00
scalar b3        = 7.5962938329E-01  
scalar seb3      = 1.9566123451E-01
scalar b4        = 1.2792483859E+00  
scalar seb4      = 6.8761936385E-01
 
qui input double(y x)
      16.08E0     1.0E0
      33.83E0     2.0E0
      65.80E0     3.0E0
      97.20E0     4.0E0
     191.55E0     5.0E0
     326.20E0     6.0E0
     386.87E0     7.0E0
     520.53E0     8.0E0
     590.03E0     9.0E0
     651.92E0    10.0E0
     724.93E0    11.0E0
     699.56E0    12.0E0
     689.96E0    13.0E0
     637.56E0    14.0E0
     717.41E0    15.0E0
end

/* The following starting values led to convergence problems:

nl ( y = {b1} / ( ( 1 + exp( {b2} - {b3}*x ) )^(1/{b4})) ), ///
	init(b1 100 b2 10 b3 1 b4 1)

*/

nl ( y = {b1} / ( ( 1 + exp( {b2} - {b3}*x ) )^(1/{b4})) ), ///
	init(b1 700 b2 5 b3 0.75 b4 1.3) eps(1e-10)

assert N    == e(N)
assert df_r == e(df_r)
assert df_m == e(df_m)

lrecomp [b1]_b[_cons] b1 [b2]_b[_cons] b2 [b3]_b[_cons] b3 [b4]_b[_cons] b4 () /*
*/ [b1]_se[_cons] seb1 [b2]_se[_cons] seb2 [b3]_se[_cons] seb3 [b4]_se[_cons] seb4 () /*
*/ e(rmse) rmse e(rss) rss

