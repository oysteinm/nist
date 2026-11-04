/* NIST/ITL StRD benchmark 

Linear Regression

Difficulty=Higher  Polynomial  k=6  N=21  Generated

Dataset Name:  Wampler-3 (wampler3.dat)

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
     Parameter          Estimate             of Estimate

        B0        1.00000000000000         2152.32624678170
        B1        1.00000000000000         2363.55173469681
        B2        1.00000000000000         779.343524331583
        B3        1.00000000000000         101.475507550350
        B4        1.00000000000000         5.64566512170752
        B5        1.00000000000000         0.112324854679312

     Residual
     Standard Deviation   2360.14502379268

     R-Squared            0.999995559025820


               Certified Analysis of Variance Table

Source of Degrees of     Sums of              Mean
Variation  Freedom       Squares             Squares           F Statistic

Regression    5      18814317208116.7    3762863441623.33   675524.458240122
Residual     15      83554268.0000000    5570284.53333333
*/

clear

scalar N        = 21
scalar df_r     = 15
scalar df_m     = 5

scalar rmse     = 2360.14502379268
scalar r2       = 0.999995559025820
scalar mss      = 18814317208116.7
scalar F        = 675524.458240122
scalar rss      = 83554268.0000000

scalar b_cons   = 1
scalar se_cons  = 2152.32624678170
scalar bx1      = 1
scalar sex1     = 2363.55173469681
scalar bx2      = 1
scalar sex2     = 779.343524331583
scalar bx3      = 1
scalar sex3     = 101.475507550350
scalar bx4      = 1
scalar sex4     = 5.64566512170752
scalar bx5      = 1
scalar sex5     = 0.112324854679312

qui input long y byte x1
              760     0
            -2042     1
             2111     2
            -1684     3
             3888     4
             1858     5
            11379     6
            17560     7
            39287     8
            64382     9
           113159    10
           175108    11
           273291    12
           400186    13
           581243    14
           811568    15
          1121004    16
          1506550    17
          2002767    18
          2611612    19
          3369180    20
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

