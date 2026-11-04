/* NIST/ITL StRD benchmark 

Linear Regression

Difficulty=Higher  Polynomial  k=6  N=21  Generated

Dataset Name:  Wampler-1 (wampler1.dat)

Procedure:     Linear Least Squares Regression

Reference:     Wampler, R. H. (1970).
               A Report of the Accuracy of Some Widely-Used Least
               Squares Computer Programs.
               Journal of the American Statistical Association, 65, pp. 549-565.

Data:          1 Response Variable (y)
               1 Predictor Variable (x)
               21 Observations
               Higher Level of Difficulty
               Generated Data

Model:         Polynomial Class
               6 Parameters (B0,B1,...,B5)

               y = B0 + B1*x + B2*(x**2) + B3*(x**3)+ B4*(x**4) + B5*(x**5)


               Certified Regression Statistics

                                          Standard Deviation
     Parameter        Estimate               of Estimate

        B0        1.00000000000000        0.000000000000000
        B1        1.00000000000000        0.000000000000000
        B2        1.00000000000000        0.000000000000000
        B3        1.00000000000000        0.000000000000000
        B4        1.00000000000000        0.000000000000000
        B5        1.00000000000000        0.000000000000000

     Residual
     Standard Deviation   0.000000000000000

     R-Squared            1.00000000000000


               Certified Analysis of Variance Table

Source of Degrees of     Sums of               Mean
Variation  Freedom       Squares              Squares           F Statistic

Regression    5      18814317208116.7     3762863441623.33       Infinity
Residual     15      0.000000000000000    0.000000000000000
*/

clear

scalar N        = 21
scalar df_r     = 15
scalar df_m     = 5

scalar rmse     = 0
scalar r2       = 1
scalar mss      = 18814317208116.7
scalar F        = .
scalar rss      = 0

scalar b_cons   = 1
scalar se_cons  = 0
scalar bx1      = 1
scalar sex1     = 0
scalar bx2      = 1
scalar sex2     = 0
scalar bx3      = 1
scalar sex3     = 0
scalar bx4      = 1
scalar sex4     = 0
scalar bx5      = 1
scalar sex5     = 0

qui input long y byte x1
                 1     0
                 6     1
                63     2
               364     3
              1365     4
              3906     5
              9331     6
             19608     7
             37449     8
             66430     9
            111111    10
            177156    11
            271453    12
            402234    13
            579195    14
            813616    15
           1118481    16
           1508598    17
           2000719    18
           2613660    19
           3368421    20
end

gen int  x2 = x1*x1
gen long x3 = x1*x2
gen long x4 = x1*x3
gen long x5 = x1*x4

reg y x1-x5
di "R-squared = " %20.15f e(r2)

assert N    == e(N)
assert df_r == e(df_r)
assert df_m == e(df_m)

lrecomp _b[_cons] b_cons _b[x1] bx1 _b[x2] bx2 /*
*/ _b[x3] bx3 _b[x4] bx4 _b[x5] bx5 () /*
*/ _se[_cons] se_cons _se[x1] sex1 _se[x2] sex2 /*
*/ _se[x3] sex3 _se[x4] sex4 _se[x5] sex5 () /*
*/ e(rmse) rmse e(r2) r2 e(mss) mss e(F) F e(rss) rss
