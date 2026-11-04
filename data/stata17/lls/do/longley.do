/* NIST/ITL StRD

Linear Regression

Difficulty=Higher  Multilinear  k=7  N=16  Observed

Dataset Name:  Longley (longley.dat)

Procedure:     Linear Least Squares Regression

Reference:     Longley, J. W. (1967).
               An Appraisal of Least Squares Programs for the
               Electronic Computer from the Viewpoint of the User.
               Journal of the American Statistical Association, 62, pp. 819-841.

Data:          1 Response Variable (y)
               6 Predictor Variable (x)
               16 Observations
               Higher Level of Difficulty
               Observed Data

Model:         Polynomial Class
               7 Parameters (B0,B1,...,B7)

               y = B0 + B1*x1 + B2*x2 + B3*x3 + B4*x4 + B5*x5 + B6*x6 + e


               Certified Regression Statistics

                                          Standard Deviation
     Parameter         Estimate              of Estimate

        B0       -3482258.63459582         890420.383607373
        B1        15.0618722713733         84.9149257747669
        B2       -0.358191792925910E-01    0.334910077722432E-01
        B3       -2.02022980381683         0.488399681651699
        B4       -1.03322686717359         0.214274163161675
        B5       -0.511041056535807E-01    0.226073200069370
        B6        1829.15146461355         455.478499142212

     Residual
     Standard Deviation   304.854073561965

     R-Squared            0.995479004577296


               Certified Analysis of Variance Table

Source of Degrees of     Sums of             Mean
Variation  Freedom       Squares            Squares          F Statistic

Regression    6      184172401.944494   30695400.3240823   330.285339234588
Residual      9      836424.055505915   92936.0061673238
*/

clear

scalar N        = 16
scalar df_r     = 9
scalar df_m     = 6

scalar rmse     = 304.854073561965
scalar r2       = 0.995479004577296
scalar mss      = 184172401.944494
scalar F        = 330.285339234588
scalar rss      = 836424.055505915

scalar b_cons   = -3482258.63459582
scalar se_cons  = 890420.383607373
scalar bx1      = 15.0618722713733
scalar sex1     = 84.9149257747669
scalar bx2      = -0.358191792925910E-01
scalar sex2     = 0.334910077722432E-01
scalar bx3      = -2.02022980381683
scalar sex3     = 0.488399681651699
scalar bx4      = -1.03322686717359
scalar sex4     = 0.214274163161675
scalar bx5      = -0.511041056535807E-01
scalar sex5     = 0.226073200069370
scalar bx6      = 1829.15146461355
scalar sex6     = 455.478499142212

qui input double (y x1 x2 x3 x4 x5 x6)
60323    83.0   234289   2356     1590    107608  1947
61122    88.5   259426   2325     1456    108632  1948
60171    88.2   258054   3682     1616    109773  1949
61187    89.5   284599   3351     1650    110929  1950
63221    96.2   328975   2099     3099    112075  1951
63639    98.1   346999   1932     3594    113270  1952
64989    99.0   365385   1870     3547    115094  1953
63761   100.0   363112   3578     3350    116219  1954
66019   101.2   397469   2904     3048    117388  1955
67857   104.6   419180   2822     2857    118734  1956
68169   108.4   442769   2936     2798    120445  1957
66513   110.8   444546   4681     2637    121950  1958
68655   112.6   482704   3813     2552    123366  1959
69564   114.2   502601   3931     2514    125368  1960
69331   115.7   518173   4806     2572    127852  1961
70551   116.9   554894   4007     2827    130081  1962
end

reg y x1-x6
di "R-squared = " %20.15f e(r2)

assert N    == e(N)
assert df_r == e(df_r)
assert df_m == e(df_m)

lrecomp _b[_cons] b_cons _b[x1] bx1 _b[x2] bx2 /*
*/ _b[x3] bx3 _b[x4] bx4 _b[x4] bx4 _b[x5] bx5 _b[x6] bx6 () /*
*/ _se[_cons] se_cons _se[x1] sex1 _se[x2] sex2 /*
*/ _se[x3] sex3 _se[x4] sex4 _se[x4] sex4 _se[x5] sex5 _se[x6] sex6 () /*
*/ e(rmse) rmse e(r2) r2 e(mss) mss e(F) F e(rss) rss


