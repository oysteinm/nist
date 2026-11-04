/* NIST/ITL StRD

Nonlinear Regression

Difficulty=Lower  Miscellaneous  k=2  N=6  Observed

Dataset Name:  DanielWood        (DanielWood.dat)

Procedure:     Nonlinear Least Squares Regression

Description:   These data and model are described in Daniel and Wood
               (1980), and originally published in E.S.Keeping,
               "Introduction to Statistical Inference," Van Nostrand
               Company, Princeton, NJ, 1962, p. 354.  The response
               variable is energy radieted from a carbon filament
               lamp per cm**2 per second, and the predictor variable
               is the absolute temperature of the filament in 1000
               degrees Kelvin.

Reference:     Daniel, C. and F. S. Wood (1980).
               Fitting Equations to Data, Second Edition.
               New York, NY:  John Wiley and Sons, pp. 428-431.


Data:          1 Response Variable  (y = energy)
               1 Predictor Variable (x = temperature)
               6 Observations
               Lower Level of Difficulty
               Observed Data

Model:         Miscellaneous Class
               2 Parameters (b1 and b2)

               y  = b1*x**b2  +  e



          Starting values                  Certified Values

        Start 1     Start 2           Parameter     Standard Deviation
  b1 =   1           0.7           7.6886226176E-01  1.8281973860E-02
  b2 =   5           4             3.8604055871E+00  5.1726610913E-02

Residual Sum of Squares:                    4.3173084083E-03
Residual Standard Deviation:                3.2853114039E-02
Degrees of Freedom:                                4
Number of Observations:                            6
*/

clear

scalar N         = 6
scalar df_r      = 4
scalar df_m      = 2

scalar rss       = 4.3173084083E-03
scalar rmse      = 3.2853114039E-02

scalar b1        = 7.6886226176E-01
scalar seb1      = 1.8281973860E-02
scalar b2        = 3.8604055871E+00
scalar seb2      = 5.1726610913E-02

qui input double (y x)
      2.138E0        1.309E0
      3.421E0        1.471E0
      3.597E0        1.490E0
      4.340E0        1.565E0
      4.882E0        1.611E0
      5.660E0        1.680E0
end

nl ( y  = {b1}*x^{b2} ), init(b1 1 b2 5) eps(1e-10) 

assert N    == e(N)
assert df_r == e(df_r)
assert df_m == e(df_m)

lrecomp [b1]_b[_cons] b1 [b2]_b[_cons] b2 () /*
*/ [b1]_se[_cons] seb1 [b2]_se[_cons] seb2 () /*
*/ e(rmse) rmse e(rss) rss
