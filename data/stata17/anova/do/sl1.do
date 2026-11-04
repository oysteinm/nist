/* NIST/ITL StRD

ANOVA

Difficulty=Lower  n_i=21  k=9  Generated

Dataset Name:   Simon-Lesage1   (Simon-Lesage1.dat)

Procedure:      Analysis of Variance


Reference:      Simon, Stephen D. and Lesage, James P. (1989).
                "Assessing the Accuracy of ANOVA Calculations in
                Statistical Software".
                Computational Statistics & Data Analysis, 8, pp. 325-332.


Data:           1 Factor
                9 Treatments
                21 Replicates/Cell
                189 Observations
                1 Constant Leading Digit
                Lower Level of Difficulty
                Generated Data


Model:          10 Parameters (mu,tau_1, ... , tau_9)
                y_{ij} = mu + tau_i + epsilon_{ij}


Certified Values:

Source of                  Sums of               Mean
Variation          df      Squares              Squares            F Statistic

Between Treatment   8 1.68000000000000E+00 2.10000000000000E-01 2.1000000000E+01
Within Treatment  180 1.80000000000000E+00 1.00000000000000E-02

                  Certified R-Squared 4.82758620689655E-01

                  Certified Residual
                  Standard Deviation  1.00000000000000E-01
*/

clear

scalar N         = 189
scalar df_r      = 180
scalar df_m      = 8

scalar mss       = 1.68
scalar F         = 21
scalar rss       = 1.8
scalar r2        = 4.82758620689655E-01
scalar rmse      = 0.1

qui input byte treat double resp
           1         1.4
           1         1.3
           1         1.5
           1         1.3
           1         1.5
           1         1.3
           1         1.5
           1         1.3
           1         1.5
           1         1.3
           1         1.5
           1         1.3
           1         1.5
           1         1.3
           1         1.5
           1         1.3
           1         1.5
           1         1.3
           1         1.5
           1         1.3
           1         1.5
           2         1.3
           2         1.2
           2         1.4
           2         1.2
           2         1.4
           2         1.2
           2         1.4
           2         1.2
           2         1.4
           2         1.2
           2         1.4
           2         1.2
           2         1.4
           2         1.2
           2         1.4
           2         1.2
           2         1.4
           2         1.2
           2         1.4
           2         1.2
           2         1.4
           3         1.5
           3         1.4
           3         1.6
           3         1.4
           3         1.6
           3         1.4
           3         1.6
           3         1.4
           3         1.6
           3         1.4
           3         1.6
           3         1.4
           3         1.6
           3         1.4
           3         1.6
           3         1.4
           3         1.6
           3         1.4
           3         1.6
           3         1.4
           3         1.6
           4         1.3
           4         1.2
           4         1.4
           4         1.2
           4         1.4
           4         1.2
           4         1.4
           4         1.2
           4         1.4
           4         1.2
           4         1.4
           4         1.2
           4         1.4
           4         1.2
           4         1.4
           4         1.2
           4         1.4
           4         1.2
           4         1.4
           4         1.2
           4         1.4
           5         1.5
           5         1.4
           5         1.6
           5         1.4
           5         1.6
           5         1.4
           5         1.6
           5         1.4
           5         1.6
           5         1.4
           5         1.6
           5         1.4
           5         1.6
           5         1.4
           5         1.6
           5         1.4
           5         1.6
           5         1.4
           5         1.6
           5         1.4
           5         1.6
           6         1.3
           6         1.2
           6         1.4
           6         1.2
           6         1.4
           6         1.2
           6         1.4
           6         1.2
           6         1.4
           6         1.2
           6         1.4
           6         1.2
           6         1.4
           6         1.2
           6         1.4
           6         1.2
           6         1.4
           6         1.2
           6         1.4
           6         1.2
           6         1.4
           7         1.5
           7         1.4
           7         1.6
           7         1.4
           7         1.6
           7         1.4
           7         1.6
           7         1.4
           7         1.6
           7         1.4
           7         1.6
           7         1.4
           7         1.6
           7         1.4
           7         1.6
           7         1.4
           7         1.6
           7         1.4
           7         1.6
           7         1.4
           7         1.6
           8         1.3
           8         1.2
           8         1.4
           8         1.2
           8         1.4
           8         1.2
           8         1.4
           8         1.2
           8         1.4
           8         1.2
           8         1.4
           8         1.2
           8         1.4
           8         1.2
           8         1.4
           8         1.2
           8         1.4
           8         1.2
           8         1.4
           8         1.2
           8         1.4
           9         1.5
           9         1.4
           9         1.6
           9         1.4
           9         1.6
           9         1.4
           9         1.6
           9         1.4
           9         1.6
           9         1.4
           9         1.6
           9         1.4
           9         1.6
           9         1.4
           9         1.6
           9         1.4
           9         1.6
           9         1.4
           9         1.6
           9         1.4
           9         1.6
end

anova resp treat

assert N    == e(N)
assert df_r == e(df_r)
assert df_m == e(df_m)

lrecomp e(F) F e(rmse) rmse e(r2) r2 e(mss) mss e(rss) rss
