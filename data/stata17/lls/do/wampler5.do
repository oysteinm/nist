/* NIST/ITL StRD benchmark

Linear Regression

Difficulty=Higher  Polynomial  k=6  N=21  Generated

Dataset Name:  Wampler-5 (wampler5.dat)

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

        B0        1.00000000000000         21523262.4678170
        B1        1.00000000000000         23635517.3469681
        B2        1.00000000000000         7793435.24331583
        B3        1.00000000000000         1014755.07550350
        B4        1.00000000000000         56456.6512170752
        B5        1.00000000000000         1123.24854679312

     Residual
     Standard Deviation   23601450.2379268

     R-Squared            0.224668921574940E-02


               Certified Analysis of Variance Table

Source of Degrees of    Sums of               Mean
Variation  Freedom      Squares              Squares           F Statistic

Regression    5    18814317208116.7      3762863441623.33 6.7552445824012241E-03
Residual     15    0.835542680000000E+16 557028453333333.
*/

clear

scalar N        = 21
scalar df_r     = 15
scalar df_m     = 5

scalar rmse     = 23601450.2379268
scalar r2       = 0.224668921574940E-02
scalar mss      = 18814317208116.7
scalar F        = 6.7552445824012241E-03
scalar rss      = 0.835542680000000E+16

scalar b_cons   = 1
scalar se_cons  = 21523262.4678170
scalar bx1      = 1
scalar sex1     = 23635517.3469681
scalar bx2      = 1
scalar sex2     = 7793435.24331583
scalar bx3      = 1
scalar sex3     = 1014755.07550350
scalar bx4      = 1
scalar sex4     = 56456.6512170752
scalar bx5      = 1
scalar sex5     = 1123.24854679312

qui input long y byte x1
             7590001     0
           -20479994     1
            20480063     2
           -20479636     3
            25231365     4
           -20476094     5
            20489331     6
           -20460392     7
            18417449     8
           -20413570     9
            20591111    10
           -20302844    11
            18651453    12
           -20077766    13
            21059195    14
           -19666384    15
            26348481    16
           -18971402    17
            22480719    18
           -17866340    19
            10958421    20
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

