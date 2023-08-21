\ENCODING=ISO-8859-1
\Problem name: example_3d_irreg_pyramid.lp

Minimize
 COST: 0 x0 + 0 x1 - 1 x2
Subject To
 GE1: -1 x0 + 0 x1 - 1 x2 >= -1.5
 GE1:  0 x0 - 1 x1 - 1 x2 >= -1.5
 GE1:  1 x0 + 0 x1 - 1 x2 >= -0.5
 GE1:  0 x0 + 1 x1 - 1 x2 >= -0.5
Bounds
 -1 <= x0 <= 1
 -1 <= x1 <= 1
 0 <= x2 <= 1
Integers
x0 x1
End
