\ENCODING=ISO-8859-1
\Problem name: example_phd.lp

Minimize
 COST: -0.01 X1 - 1 X3
Subject To
 GE1:  X1 + X3 <= 1
 GE2:  0.5 X2 + X3 <= 0.5
 GE3:  X1 + X2 >= 0.5
Bounds
 X1 >= 0
 X2 >= 0
Integers
X1 X2
End
