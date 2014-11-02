*NAME:         small
*ROWS:         2
*COLUMNS:      3
*INTEGER:      1
*NONZERO:      4
*BEST SOLN:    568.101 (opt)
*LP SOLN:      149.589
*SOURCE:       Aykut Bulut, aykutblt@gmail.com
*APPLICATION:  created to test Osi's simplexrelated methods.
*COMMENTS:     created to test Osi's simplexrelated methods.
*
NAME          small
ROWS
 N  cost
 E  row1
COLUMNS
    x1        cost                 1   row1                 1
    x2        cost                 1   row1                 2
RHS
    rhs1      row1                 2
BOUNDS
 LO bnd1      x1
 LO bnd1      x2
ENDATA
