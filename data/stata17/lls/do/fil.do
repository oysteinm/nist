/* NIST/ITL StRD

Linear Regression

Difficulty=Higher Polynomial k=11 N=82 Observed

Dataset Name:  Filippelli (filippelli.dat)

Procedure:     Linear Least Squares Regression

Reference:     Filippelli, A., NIST.

Data:          1 Response Variable (y)
               1 Predictor Variable (x)
               82 Observations
               Higher Level of Difficulty
               Observed Data

Model:         Polynomial Class
               11 Parameters (B0,B1,...,B10)

               y = B0 + B1*x + B2*(x**2) + ... + B9*(x**9) + B10*(x**10) + e


               Certified Regression Statistics

                                            Standard Deviation
     Parameter         Estimate                of Estimate

        B0        -1467.48961422980         298.084530995537
        B1        -2772.17959193342         559.779865474950
        B2        -2316.37108160893         466.477572127796
        B3        -1127.97394098372         227.204274477751
        B4        -354.478233703349         71.6478660875927
        B5        -75.1242017393757         15.2897178747400
        B6        -10.8753180355343         2.23691159816033
        B7        -1.06221498588947         0.221624321934227
        B8        -0.670191154593408E-01    0.142363763154724E-01
        B9        -0.246781078275479E-02    0.535617408889821E-03
        B10       -0.402962525080404E-04    0.896632837373868E-05

     Residual
     Standard Deviation   0.334801051324544E-02

     R-Squared            0.996727416185620


               Certified Analysis of Variance Table

Source of Degrees of     Sums of                 Mean
Variation  Freedom       Squares                Squares           F Statistic

Regression   10     0.242391619837339     0.242391619837339E-01 2162.43954511489
Residual     71     0.795851382172941E-03 0.112091743968020E-04
*/

clear

scalar N        = 82
scalar df_r     = 71
scalar df_m     = 10

scalar rmse     = 0.334801051324544E-02
scalar r2       = 0.996727416185620
scalar mss      = 0.242391619837339
scalar F        = 2162.43954511489
scalar rss      = 0.795851382172941E-03

scalar b_cons   = -1467.48961422980
scalar se_cons  = 298.084530995537
scalar bx1      = -2772.17959193342
scalar sex1     = 559.779865474950
scalar bx2      = -2316.37108160893
scalar sex2     = 466.477572127796
scalar bx3      = -1127.97394098372
scalar sex3     = 227.204274477751
scalar bx4      = -354.478233703349
scalar sex4     = 71.6478660875927
scalar bx5      = -75.1242017393757
scalar sex5     = 15.2897178747400
scalar bx6      = -10.8753180355343
scalar sex6     = 2.23691159816033
scalar bx7      = -1.06221498588947
scalar sex7     = 0.221624321934227
scalar bx8      = -0.670191154593408E-01
scalar sex8     = 0.142363763154724E-01
scalar bx9      = -0.246781078275479E-02
scalar sex9     = 0.535617408889821E-03
scalar bx10     = -0.402962525080404E-04
scalar sex10    = 0.896632837373868E-05

qui input double (y x1)
            0.8116   -6.860120914
            0.9072   -4.324130045
            0.9052   -4.358625055
            0.9039   -4.358426747
            0.8053   -6.955852379
            0.8377   -6.661145254
            0.8667   -6.355462942
            0.8809   -6.118102026
            0.7975   -7.115148017
            0.8162   -6.815308569
            0.8515   -6.519993057
            0.8766   -6.204119983
            0.8885   -5.853871964
            0.8859   -6.109523091
            0.8959   -5.79832982
            0.8913   -5.482672118
            0.8959   -5.171791386
            0.8971   -4.851705903
            0.9021   -4.517126416
            0.909    -4.143573228
            0.9139   -3.709075441
            0.9199   -3.499489089
            0.8692   -6.300769497
            0.8872   -5.953504836
            0.89     -5.642065153
            0.891    -5.031376979
            0.8977   -4.680685696
            0.9035   -4.329846955
            0.9078   -3.928486195
            0.7675   -8.56735134
            0.7705   -8.363211311
            0.7713   -8.107682739
            0.7736   -7.823908741
            0.7775   -7.522878745
            0.7841   -7.218819279
            0.7971   -6.920818754
            0.8329   -6.628932138
            0.8641   -6.323946875
            0.8804   -5.991399828
            0.7668   -8.781464495
            0.7633   -8.663140179
            0.7678   -8.473531488
            0.7697   -8.247337057
            0.77     -7.971428747
            0.7749   -7.676129393
            0.7796   -7.352812702
            0.7897   -7.072065318
            0.8131   -6.774174009
            0.8498   -6.478861916
            0.8741   -6.159517513
            0.8061   -6.835647144
            0.846    -6.53165267
            0.8751   -6.224098421
            0.8856   -5.910094889
            0.8919   -5.598599459
            0.8934   -5.290645224
            0.894    -4.974284616
            0.8957   -4.64454848
            0.9047   -4.290560426
            0.9129   -3.885055584
            0.9209   -3.408378962
            0.9219   -3.13200249
            0.7739   -8.726767166
            0.7681   -8.66695597
            0.7665   -8.511026475
            0.7703   -8.165388579
            0.7702   -7.886056648
            0.7761   -7.588043762
            0.7809   -7.283412422
            0.7961   -6.995678626
            0.8253   -6.691862621
            0.8602   -6.392544977
            0.8809   -6.067374056
            0.8301   -6.684029655
            0.8664   -6.378719832
            0.8834   -6.065855188
            0.8898   -5.752272167
            0.8964   -5.132414673
            0.8963   -4.811352704
            0.9074   -4.098269308
            0.9119   -3.66174277
            0.9228   -3.2644011
end

gen double x2 = x1*x1
gen double x3 = x1*x2
gen double x4 = x1*x3
gen double x5 = x1*x4
gen double x6 = x1*x5
gen double x7 = x1*x6
gen double x8 = x1*x7
gen double x9 = x1*x8
gen double x10 = x1*x9

reg y x1-x10
di "R-squared = " %20.15f e(r2)

assert N    == e(N)
cap noi assert df_r == e(df_r)
cap noi assert df_m == e(df_m)

lrecomp _b[_cons] b_cons _b[x1] bx1 _b[x2] bx2 /*
*/ _b[x3] bx3 _b[x4] bx4 _b[x5] bx5 _b[x6] bx6 /*
*/ _b[x7] bx7 _b[x8] bx8 _b[x9] bx9 _b[x10] bx10 () /*
*/ _se[_cons] se_cons _se[x1] sex1 _se[x2] sex2 /*
*/ _se[x3] sex3 _se[x4] sex4 _se[x5] sex5 _se[x6] sex6 /*
*/ _se[x7] sex7 _se[x8] sex8 _se[x9] sex9 _se[x10] sex10 () /*
*/ e(rmse) rmse e(r2) r2 e(mss) mss e(F) F e(rss) rss

* Do alternative calculation of rmse

predict double res, residual
gen double ss = res*res
qui summarize ss
scalar rmsea = sqrt(r(sum)/e(df_r))
di _n in gr "RMSE standard       = " in ye %22.15e rmse _n /*
*/ in gr    "RMSE from residual  = " in ye %22.15e rmsea _n /*
*/ in gr    "abs(std. - res.)    = " in ye %22.15e abs(rmsea-rmse) _n /*
*/ in gr    "RMSE from regress   = " in ye %22.15e e(rmse) _n /*
*/ in gr    "abs(res. - regress) = " in ye %22.15e abs(rmsea-e(rmse))

* Use -orthog-

orthog x*, gen(u*) mat(R)
reg y u*
matrix b = e(b)*inv(R)'
lrecomp b[1,11] b_cons b[1,1] bx1 b[1,2] bx2 /*
*/ b[1,3] bx3 b[1,4] bx4 b[1,5] bx5 b[1,6] bx6 /*
*/ b[1,7] bx7 b[1,8] bx8 b[1,9] bx9 b[1,10] bx10 () /*
*/ e(rmse) rmse e(r2) r2 e(mss) mss e(F) F e(rss) rss

