*NAME:         small
*ROWS:         4
*COLUMNS:      3
*INTEGER:      1
*NONZERO:      7
*BEST SOLN:
*LP SOLN:
*SOURCE:       Aykut Bulut, aykutblt@gmail.com
*APPLICATION:  created to test Osi's simplexrelated methods.
*COMMENTS:     created to test Osi's simplexrelated methods.
*
NAME          small
ROWS
 N  cost
 G  row1
 L  row2
 E  row3
COLUMNS
    xone      cost                 2   row1                 1
    xone      row2                -1   row3                 1
    xtwo      cost                 1   row1                 1
    xtwo      row3                 -1
    xthree    cost                -1   row2                 1
    xthree    row3                 1
RHS
    rhs1      row1                 2   row2                 2
    rhs1      row3                 2
BOUNDS
 LI bnd1      xone
 LO bnd1      xtwo
 LO bnd1      xthree
ENDATA
