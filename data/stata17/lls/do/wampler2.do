/* NIST/ITL StRD benchmark 

Linear Regression

Difficulty=Higher  Polynomial  k=6  N=21  Generated

Dataset Name:  Wampler-2 (wampler2.dat)

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
     Parameter         Estimate               of Estimate

        B0        1.00000000000000         0.000000000000000
        B1        0.100000000000000        0.000000000000000
        B2        0.100000000000000E-01    0.000000000000000
        B3        0.100000000000000E-02    0.000000000000000
        B4        0.100000000000000E-03    0.000000000000000
        B5        0.100000000000000E-04    0.000000000000000

     Residual
     Standard Deviation   0.000000000000000
     R-Squared            1.00000000000000


               Certified Analysis of Variance Table

Source of Degrees of     Sums of               Mean
Variation  Freedom       Squares              Squares           F Statistic

Regression   5       6602.91858365167     1320.58371673033       Infinity
Residual    15       0.000000000000000    0.000000000000000
*/

clear

scalar N        = 21
scalar df_r     = 15
scalar df_m     = 5

scalar rmse     = 0
scalar r2       = 1
scalar mss      = 6602.91858365167
scalar F        = .
scalar rss      = 0

scalar b_cons   = 1
scalar se_cons  = 0
scalar bx1      = 1e-1
scalar sex1     = 0
scalar bx2      = 1e-2
scalar sex2     = 0
scalar bx3      = 1e-3
scalar sex3     = 0
scalar bx4      = 1e-4
scalar sex4     = 0
scalar bx5      = 1e-5
scalar sex5     = 0

qui input double y byte x1
            1.00000    0
            1.11111    1
            1.24992    2
            1.42753    3
            1.65984    4
            1.96875    5
            2.38336    6
            2.94117    7
            3.68928    8
            4.68559    9
            6.00000   10
            7.71561   11
            9.92992   12
           12.75603   13
           16.32384   14
           20.78125   15
           26.29536   16
           33.05367   17
           41.26528   18
           51.16209   19
           63.00000   20
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

