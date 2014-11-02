*NAME:         mir_example
*ROWS:
*COLUMNS:
*INTEGER:
*NONZERO:
*BEST SOLN:
*LP SOLN:
*SOURCE:       Aykut Bulut, aykutblt@gmail.com
*APPLICATION:  created to test Atamturk and Narayanan's conic MIR cuts.
*COMMENTS:     created to test Atamturk and Narayanan's conic MIR cuts.
*
* min 5x1+x2+2x3+x4
* s.t.
* x1+x2    +x4 >= 3.5
* x1   -x3 +x4 <= 6.5
* x1+x2+x3     =  4.5
* (x1,x2,x3) in L^3
* x2 integer
* x1,x3,x4 >= 0
*
*
NAME          mir_example
ROWS
 N  cost
 G  row1
 L  row2
 E  row3
COLUMNS
    x1        cost                 5   row1                 1
    x1        row2                 1   row3                 1
    x2        cost                 1   row1                 1
    x2        row3                 1
    x3        cost                 2   row2                -1
    x3        row3                 1
    x4        cost                 1   row1                -1
    x4        row2                 1
RHS
    rhs1      row1                 3.5 row2                 6.5
    rhs1      row3                 4.5
BOUNDS
 LO bnd1      x1
 LI bnd1      x2
 LO bnd1      x3
 LO bnd1      x4
CSECTION      k1        0              QUAD
    x1
    x2
    x3
ENDATA
