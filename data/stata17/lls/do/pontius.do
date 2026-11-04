/* NIST/ITL StRD

Linear Regression

Difficulty=Lower  Quadratic  k=3  N=40  Observed

Dataset Name:  Pontius

Procedure:     Linear Least Squares Regression

Reference:     Pontius, P., NIST.
               Load Cell Calibration.

Data:          1 Response Variable (y)
               1 Predictor Variable (x)
               40 Observations
               Lower Level of Difficulty
               Observed Data

Model:         Quadratic Class
               3 Parameters (B0,B1,B2)
               y = B0 + B1*x + B2*(x**2)


               Certified Regression Statistics

                                          Standard Deviation
     Parameter          Estimate             of Estimate

        B0      0.673565789473684E-03    0.107938612033077E-03
        B1      0.732059160401003E-06    0.157817399981659E-09
        B2     -0.316081871345029E-14    0.486652849992036E-16

     Residual
     Standard Deviation   0.205177424076185E-03

     R-Squared            0.999999900178537


               Certified Analysis of Variance Table

Source of Degrees of     Sums of               Mean
Variation  Freedom       Squares              Squares             F Statistic

Regression    2     15.6040343244198      7.80201716220991      185330865.995752
Residual     37     0.155761768796992E-05 0.420977753505385E-07
*/

clear

scalar N        = 40
scalar df_r     = 37
scalar df_m     = 2

scalar rmse     = 0.205177424076185E-03
scalar r2       = 0.999999900178537
scalar mss      = 15.6040343244198
scalar F        = 185330865.995752
scalar rss      = 0.155761768796992E-05

scalar b_cons   = 0.673565789473684E-03
scalar se_cons  = 0.107938612033077E-03
scalar bx1      = 0.732059160401003E-06
scalar sex1     = 0.157817399981659E-09
scalar bx2      = -0.316081871345029E-14
scalar sex2     = 0.486652849992036E-16

qui input double (y x1)
         .11019        150000
         .21956        300000
         .32949        450000
         .43899        600000
         .54803        750000
         .65694        900000
         .76562       1050000
         .87487       1200000
         .98292       1350000
        1.09146       1500000
        1.20001       1650000
        1.30822       1800000
        1.41599       1950000
        1.52399       2100000
        1.63194       2250000
        1.73947       2400000
        1.84646       2550000
        1.95392       2700000
        2.06128       2850000
        2.16844       3000000
         .11052        150000
         .22018        300000
         .32939        450000
         .43886        600000
         .54798        750000
         .65739        900000
         .76596       1050000
         .87474       1200000
         .98300       1350000
        1.09150       1500000
        1.20004       1650000
        1.30818       1800000
        1.41613       1950000
        1.52408       2100000
        1.63159       2250000
        1.73965       2400000
        1.84696       2550000
        1.95445       2700000
        2.06177       2850000
        2.16829       3000000
end

gen double x2 = x1*x1

reg y x1 x2
di "R-squared = " %20.15f e(r2)

assert N    == e(N)
assert df_r == e(df_r)
assert df_m == e(df_m)

lrecomp _b[_cons] b_cons _b[x1] bx1 _b[x2] bx2 () /*
*/ _se[_cons] se_cons _se[x1] sex1 _se[x2] sex2 () /*
*/ e(rmse) rmse e(r2) r2 e(mss) mss e(F) F e(rss) rss

