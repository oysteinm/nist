/* NIST/ITL StRD

Univariate Summary Statistics

Difficulty=Lower N=100 Observed

File Name:     michelso.dat

Dataset Name:  Speed of Light (in millions of meters per second)

Description:   This is an observed/"real world" data set
               consisting of 100 measurements of the
               speed of light in air.  This classic experiment
               was carried out by Michelson is 1879.
               We here use this data to test accuracy
               in summary statistics calculations.

Stat Category: Univariate: Summary Statistics

Reference:     Dorsey, Ernest N. (1944). The Velocity of Light.
               Transactions of the American Philiosophical
               Society, Volume 34, Part 1, Pages 1-110, Table 22.

               y    = mu + e

                                                  Certified Values
Sample Mean                                ybar:  299.852400000000
Sample Standard Deviation (denom. = n-1)      s:  0.0790105478190518
Sample Autocorrelation Coefficient (lag 1) r(1):  0.535199668621283

Number of Observations:                             100
*/

clear

scalar N    = 100
scalar mean = 299.852400000000
scalar sd   = 0.0790105478190518
scalar rho1 = 0.535199668621283

qui input double y
  299.85
  299.74
  299.90
  300.07
  299.93
  299.85
  299.95
  299.98
  299.98
  299.88
  300.00
  299.98
  299.93
  299.65
  299.76
  299.81
  300.00
  300.00
  299.96
  299.96
  299.96
  299.94
  299.96
  299.94
  299.88
  299.80
  299.85
  299.88
  299.90
  299.84
  299.83
  299.79
  299.81
  299.88
  299.88
  299.83
  299.80
  299.79
  299.76
  299.80
  299.88
  299.88
  299.88
  299.86
  299.72
  299.72
  299.62
  299.86
  299.97
  299.95
  299.88
  299.91
  299.85
  299.87
  299.84
  299.84
  299.85
  299.84
  299.84
  299.84
  299.89
  299.81
  299.81
  299.82
  299.80
  299.77
  299.76
  299.74
  299.75
  299.76
  299.91
  299.92
  299.89
  299.86
  299.88
  299.72
  299.84
  299.85
  299.85
  299.78
  299.89
  299.84
  299.78
  299.81
  299.76
  299.81
  299.79
  299.81
  299.82
  299.85
  299.87
  299.87
  299.81
  299.74
  299.81
  299.94
  299.95
  299.80
  299.81
  299.87
end

summarize y
assert scalar(N) == r(N)

lrecomp r(mean) mean sqrt(r(Var)) sd

gen time = _n
tsset time
corrgram y, lag(1)

lrecomp r(ac1) rho1
