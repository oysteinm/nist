/* NIST/ITL StRD

Linear Regression

Difficulty=Average  Linear k=1  N=3  Generated

Dataset Name:  Line Through Origin-2 (nointercept2.dat)

Procedure:     Linear Least Squares Regression

Reference:     Eberhardt, K., NIST.

Data:          1 Response Variable (y)
               1 Predictor Variable (x)
               3 Observations
               Average Level of Difficulty
               Generated Data

Model:         Linear Class
               1 Parameter (B1)

               y = B1*x + e


               Certified Regression Statistics

                                          Standard Deviation
     Parameter          Estimate             of Estimate

        B1          0.727272727272727    0.420827318078432E-01

     Residual
     Standard Deviation   0.369274472937998

     R-Squared            0.993348115299335


               Certified Analysis of Variance Table

Source of Degrees of    Sums of             Mean
Variation  Freedom      Squares            Squares           F Statistic

Regression   1     40.7272727272727    40.7272727272727   298.6666666666667
Residual     2     0.272727272727273   0.136363636363636
*/

clear

scalar N        = 3
scalar df_r     = 2
scalar df_m     = 1

scalar rmse     = 0.369274472937998
scalar r2       = 0.993348115299335
scalar mss      = 40.7272727272727
scalar F        = 298.6666666666667
scalar rss      = 0.272727272727273

scalar bx       = 0.727272727272727
scalar sex      = 0.420827318078432E-01

qui input int (y x)
         3       4
         4       5
         4       6
end

reg y x, nocons
di "R-squared = " %20.15f e(r2)

assert N    == e(N)
assert df_r == e(df_r)
assert df_m == e(df_m)

lrecomp _b[x] bx () _se[x] sex () /*
*/ e(rmse) rmse e(r2) r2 e(mss) mss e(F) F e(rss) rss
