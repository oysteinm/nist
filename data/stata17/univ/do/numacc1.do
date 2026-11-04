/* NIST/ITL StRD

Univariate Summary Statistics

Difficulty=Lower N=3 Generated

File Name:     numacc1.dat

Dataset Name:  Numerical-Accuracy-1

Description:   This is a constructed/fabricated data set
               to test accuracy in summary statistic calculations.
               The numbers are large (8-digit integers) and
               differ only in the last decimal place.
               Note--by construction, this data set has
                     sample mean                  = 10000002 (exact)
                     sample standard deviation    =        1 (exact)
                     sample autocorrelation coef. =     -0.5 (exact)

Stat Category: Univariate: Summary Statistics

Reference:     Simon, Stephen D. and Lesage, James P. (1989).
               Assessing the Accuracy of ANOVA Caluclations
               in Statistical Software", Computational
               Statistics & data Analysis, 8, pp. 325-332.

Data:          Constructed
               1    Response           : y
               0    Predictors
               3    Observations

Model:         Lower Level of Difficulty
               2    Parameters         : mu, sigma
               1    Response Variable  : y
               0    Predictor Variables

               y    = mu + e



                                                  Certified Values
Sample Mean                                ybar:  10000002 (exact)
Sample Standard Deviation (denom. = n-1)      s:         1 (exact)
Sample Autocorrelation Coefficient (lag 1) r(1):      -0.5 (exact)

Number of Observations:                                       3
*/

clear

scalar N    = 3
scalar mean = 10000002
scalar sd   = 1
scalar rho1 = -0.5

qui input long y
10000001
10000003
10000002
end

summarize y
assert scalar(N) == r(N)

lrecomp r(mean) mean sqrt(r(Var)) sd

gen time = _n
tsset time
corrgram y, lag(1)

lrecomp r(ac1) rho1
