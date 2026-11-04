/* NIST/ITL StRD benchmark

Linear Regression

Difficulty=Higher  Polynomial  k=6  N=21  Generated

Dataset Name:  Wampler-4 (wampler4.dat)

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

        B0        1.00000000000000         215232.624678170
        B1        1.00000000000000         236355.173469681
        B2        1.00000000000000         77934.3524331583
        B3        1.00000000000000         10147.5507550350
        B4        1.00000000000000         564.566512170752
        B5        1.00000000000000         11.2324854679312

     Residual
     Standard Deviation   236014.502379268

     R-Squared            0.957478440825662


               Certified Analysis of Variance Table

Source of Degrees of     Sums of                Mean
Variation  Freedom       Squares               Squares          F Statistic

Regression    5      18814317208116.7      3762863441623.33   67.5524458240122
Residual     15      835542680000.000      55702845333.3333
*/

clear

scalar N        = 21
scalar df_r     = 15
scalar df_m     = 5

scalar rmse     = 236014.502379268
scalar r2       = 0.957478440825662
scalar mss      = 18814317208116.7
scalar F        = 67.5524458240122
scalar rss      = 835542680000.000

scalar b_cons   = 1
scalar se_cons  = 215232.624678170
scalar bx1      = 1
scalar sex1     = 236355.173469681
scalar bx2      = 1
scalar sex2     = 77934.3524331583
scalar bx3      = 1
scalar sex3     = 10147.5507550350
scalar bx4      = 1
scalar sex4     = 564.566512170752
scalar bx5      = 1
scalar sex5     = 11.2324854679312

qui input long y byte x1
              75901    0
            -204794    1
             204863    2
            -204436    3
             253665    4
            -200894    5
             214131    6
            -185192    7
             221249    8
            -138370    9
             315911   10
             -27644   11
             455253   12
             197434   13
             783995   14
             608816   15
            1370781   16
            1303798   17
            2205519   18
            2408860   19
            3444321   20
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
