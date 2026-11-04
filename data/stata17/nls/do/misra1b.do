/* NIST/ITL StRD

Nonlinear Regression

Difficulty=Lower  Miscellaneous  k=2  N=14  Observed

Dataset Name:  Misra1b           (Misra1b.dat)

Procedure:     Nonlinear Least Squares Regression

Description:   These data are the result of a NIST study regarding
               dental research in monomolecular adsorption.  The
               response variable is volume, and the predictor
               variable is pressure.

Reference:     Misra, D., NIST (1978).
               Dental Research Monomolecular Adsorption Study.


Data:          1 Response  (y = volume)
               1 Predictor (x = pressure)
               14 Observations
               Lower Level of Difficulty
               Observed Data

Model:         Miscellaneous Class
               2 Parameters (b1 and b2)

               y = b1 * (1-(1+b2*x/2)**(-2))  +  e


          Starting values                  Certified Values

        Start 1     Start 2           Parameter     Standard Deviation
  b1 =   500         300           3.3799746163E+02  3.1643950207E+00
  b2 =     0.0001      0.0002      3.9039091287E-04  4.2547321834E-06

Residual Sum of Squares:                    7.5464681533E-02
Residual Standard Deviation:                7.9301471998E-02
Degrees of Freedom:                                12
Number of Observations:                            14
*/

clear

scalar N         = 14
scalar df_r      = 12
scalar df_m      = 2

scalar rss       = 7.5464681533E-02
scalar rmse      = 7.9301471998E-02

scalar b1        = 3.3799746163E+02
scalar seb1      = 3.1643950207E+00
scalar b2        = 3.9039091287E-04
scalar seb2      = 4.2547321834E-06

qui input double (y x)
      10.07E0      77.6E0
      14.73E0     114.9E0
      17.94E0     141.1E0
      23.93E0     190.8E0
      29.61E0     239.9E0
      35.18E0     289.0E0
      40.02E0     332.8E0
      44.82E0     378.4E0
      50.76E0     434.8E0
      55.05E0     477.3E0
      61.01E0     536.8E0
      66.40E0     593.1E0
      75.47E0     689.1E0
      81.78E0     760.0E0
end

nl ( y = {b1} * (1 - (1 + {b2}*x/2)^(-2)) ), ///
	init(b1 500 b2 0.0001) eps(1e-10)

assert N    == e(N)
assert df_r == e(df_r)
assert df_m == e(df_m)

lrecomp [b1]_b[_cons] b1 [b2]_b[_cons] b2 () /*
*/ [b1]_se[_cons] seb1 [b2]_se[_cons] seb2 () /*
*/ e(rmse) rmse e(rss) rss
