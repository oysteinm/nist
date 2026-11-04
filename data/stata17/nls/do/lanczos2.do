/* NIST/ITL StRD
Dataset Name:  Lanczos2          (Lanczos2.dat)

File Format:   ASCII
               Starting Values   (lines 41 to 46)
               Certified Values  (lines 41 to 51)
               Data              (lines 61 to 84)

Procedure:     Nonlinear Least Squares Regression

Description:   These data are taken from an example discussed in
               Lanczos (1956).  The data were generated to 6-digits
               of accuracy using
               f(x) = 0.0951*exp(-x) + 0.8607*exp(-3*x) 
                                     + 1.5576*exp(-5*x).


Reference:     Lanczos, C. (1956).
               Applied Analysis.
               Englewood Cliffs, NJ:  Prentice Hall, pp. 272-280.




Data:          1 Response  (y)
               1 Predictor (x)
               24 Observations
               Average Level of Difficulty
               Generated Data
 
Model:         Exponential Class
               6 Parameters (b1 to b6)
 
               y = b1*exp(-b2*x) + b3*exp(-b4*x) + b5*exp(-b6*x)  +  e
 


          Starting values                  Certified Values

        Start 1     Start 2           Parameter     Standard Deviation
  b1 =   1.2         0.5           9.6251029939E-02  6.6770575477E-04
  b2 =   0.3         0.7           1.0057332849E+00  3.3989646176E-03
  b3 =   5.6         3.6           8.6424689056E-01  1.7185846685E-03
  b4 =   5.5         4.2           3.0078283915E+00  4.1707005856E-03
  b5 =   6.5         4             1.5529016879E+00  2.3744381417E-03
  b6 =   7.6         6.3           5.0028798100E+00  1.3958787284E-03

Residual Sum of Squares:                    2.2299428125E-11
Residual Standard Deviation:                1.1130395851E-06
Degrees of Freedom:                                18
Number of Observations:                            24
*/

clear

scalar N         = 24
scalar df_r      = 18
scalar df_m      = 6

scalar rss       = 2.2299428125E-11
scalar rmse      = 1.1130395851E-06


scalar b1        = 9.6251029939E-02  
scalar seb1      = 6.6770575477E-04
scalar b2        = 1.0057332849E+00  
scalar seb2      = 3.3989646176E-03
scalar b3        = 8.6424689056E-01  
scalar seb3      = 1.7185846685E-03
scalar b4        = 3.0078283915E+00  
scalar seb4      = 4.1707005856E-03
scalar b5        = 1.5529016879E+00  
scalar seb5      = 2.3744381417E-03
scalar b6        = 5.0028798100E+00  
scalar seb6      = 1.3958787284E-03

qui input double(y x)
       2.51340E+00  0.00000E+00
       2.04433E+00  5.00000E-02
       1.66840E+00  1.00000E-01
       1.36642E+00  1.50000E-01
       1.12323E+00  2.00000E-01
       9.26890E-01  2.50000E-01
       7.67934E-01  3.00000E-01
       6.38878E-01  3.50000E-01
       5.33784E-01  4.00000E-01
       4.47936E-01  4.50000E-01
       3.77585E-01  5.00000E-01
       3.19739E-01  5.50000E-01
       2.72013E-01  6.00000E-01
       2.32497E-01  6.50000E-01
       1.99659E-01  7.00000E-01
       1.72270E-01  7.50000E-01
       1.49341E-01  8.00000E-01
       1.30070E-01  8.50000E-01
       1.13812E-01  9.00000E-01
       1.00042E-01  9.50000E-01
       8.83321E-02  1.00000E+00
       7.83354E-02  1.05000E+00
       6.97669E-02  1.10000E+00
       6.23931E-02  1.15000E+00
end

nl ( y = {b1}*exp(-{b2}*x) + {b3}*exp(-{b4}*x) + {b5}*exp(-{b6}*x) ), ///
	init(b1 1.2 b2 0.3 b3 5.6 b4 5.5 b5 6.5 b6 7.6) eps(1e-10)

assert N    == e(N)
assert df_r == e(df_r)
assert df_m == e(df_m)

lrecomp [b1]_b[_cons] b1 [b2]_b[_cons] b2 [b3]_b[_cons] b3 /*
*/ [b4]_b[_cons] b4 [b5]_b[_cons] b5 [b6]_b[_cons] b6 () /*
*/ [b1]_se[_cons] seb1 [b2]_se[_cons] seb2 [b3]_se[_cons] seb3 /*
*/ [b4]_se[_cons] seb4 [b5]_se[_cons] seb5 [b6]_se[_cons] seb6 () /*
*/ e(rmse) rmse e(rss) rss


