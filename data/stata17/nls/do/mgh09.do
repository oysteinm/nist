/* NIST/ITL StRD
Dataset Name:  MGH09             (MGH09.dat)

File Format:   ASCII
               Starting Values   (lines 41 to 44)
               Certified Values  (lines 41 to 49)
               Data              (lines 61 to 71)

Procedure:     Nonlinear Least Squares Regression

Description:   This problem was found to be difficult for some very 
               good algorithms.  There is a local minimum at (+inf,
               -14.07..., -inf, -inf) with final sum of squares 
               0.00102734....

               See More, J. J., Garbow, B. S., and Hillstrom, K. E. 
               (1981).  Testing unconstrained optimization software.
               ACM Transactions on Mathematical Software. 7(1): 
               pp. 17-41.

Reference:     Kowalik, J.S., and M. R. Osborne, (1978).  
               Methods for Unconstrained Optimization Problems.  
               New York, NY:  Elsevier North-Holland.

Data:          1 Response  (y)
               1 Predictor (x)
               11 Observations
               Higher Level of Difficulty
               Generated Data
 
Model:         Rational Class (linear/quadratic)
               4 Parameters (b1 to b4)
 
               y = b1*(x**2+x*b2) / (x**2+x*b3+b4)  +  e
 

 
          Starting values                  Certified Values

        Start 1     Start 2           Parameter     Standard Deviation
  b1 =   25          0.25          1.9280693458E-01  1.1435312227E-02
  b2 =   39          0.39          1.9128232873E-01  1.9633220911E-01
  b3 =   41.5        0.415         1.2305650693E-01  8.0842031232E-02
  b4 =   39          0.39          1.3606233068E-01  9.0025542308E-02

Residual Sum of Squares:                    3.0750560385E-04
Residual Standard Deviation:                6.6279236551E-03
Degrees of Freedom:                                7
Number of Observations:                           11
*/

clear

scalar N         = 11
scalar df_r      = 7
scalar df_m      = 4

scalar rss       = 3.0750560385E-04
scalar rmse      = 6.6279236551E-03

scalar b1        = 1.9280693458E-01  
scalar seb1      = 1.1435312227E-02
scalar b2        = 1.9128232873E-01  
scalar seb2      = 1.9633220911E-01
scalar b3        = 1.2305650693E-01  
scalar seb3      = 8.0842031232E-02
scalar b4        = 1.3606233068E-01  
scalar seb4      = 9.0025542308E-02

qui input double(y x)
       1.957000E-01    4.000000E+00
       1.947000E-01    2.000000E+00
       1.735000E-01    1.000000E+00
       1.600000E-01    5.000000E-01
       8.440000E-02    2.500000E-01
       6.270000E-02    1.670000E-01
       4.560000E-02    1.250000E-01
       3.420000E-02    1.000000E-01
       3.230000E-02    8.330000E-02
       2.350000E-02    7.140000E-02
       2.460000E-02    6.250000E-02
end

/* The following starting values led to convergence problems:

nl ( y = {b1}*(x^2 + x*{b2}) / (x^2 + x*{b3} + {b4}) ), ///
	init(b1 25 b2 39 b3 41.5 b4 39)

*/

nl ( y = {b1}*(x^2 + x*{b2}) / (x^2 + x*{b3} + {b4}) ), ///
	init(b1 0.25 b2 0.39 b3 0.415 b4 0.39) eps(1e-10)

assert N    == e(N)
assert df_r == e(df_r)
assert df_m == e(df_m)

lrecomp [b1]_b[_cons] b1 [b2]_b[_cons] b2 [b3]_b[_cons] b3 [b4]_b[_cons] b4 () /*
*/ [b1]_se[_cons] seb1 [b2]_se[_cons] seb2 [b3]_se[_cons] seb3 [b4]_se[_cons] seb4 () /*
*/ e(rmse) rmse e(rss) rss

