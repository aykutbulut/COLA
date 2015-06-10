#include <OsiConicSolverInterface.hpp>
#include <ColaModel.hpp>
#include <CoinPackedMatrix.hpp>
#include <CoinPackedVector.hpp>

// the problem is
// min x1 + x2
// s.t.
// || 2x1 - x2       + 2x4  - 1 || <= [-1 -1 0 1] x1 - 2
// || x1 + x2 + 3x3           2 ||                x2
//                                                x3
//                                                x4
//
//

void fill_A(CoinPackedMatrix *& mat);
void generate_conic_constraint(CoinPackedMatrix *& A, CoinPackedVector *& b,
                               CoinPackedVector *& d, double & h);
void generate_linear_constraints(CoinPackedMatrix *& matrix, double *& col_lb,
                                 double *& col_ub, double *& obj,
                                 double *& row_lb, double *& row_ub);

// linear constraint set is empty for now.
void generate_linear_constraints(CoinPackedMatrix *& matrix, double *& col_lb,
                                 double *& col_ub, double *& obj,
                                 double *& row_lb, double *& row_ub) {
  matrix = new CoinPackedMatrix(false, 0, 0);
  // int row1_ind[] = {0,1};
  // double row1_val[] = {1.0, 1.0};
  // matrix->appendRow(CoinPackedVector(2, row1_ind, row1_val));
  // // =====row2
  // int row2_ind[] = {1};
  // double row2_val[] = {2.0};
  // matrix->appendRow(CoinPackedVector(1, row2_ind, row2_val));
  // set column lower bounds
  col_lb = NULL;
  col_ub = NULL;
  obj = new double[4]();
  obj[0] = 1.0;
  obj[2] = 1.0;
  row_lb = NULL;
  row_ub = NULL;
}

void fill_A(CoinPackedMatrix *& A) {
  A = new CoinPackedMatrix(false, 0, 0);
  int index1[] = {0, 1, 3};
  int index2[] = {0, 1, 2};
  double val1[] = {2.0, -1.0, 2.0};
  double val2[] = {1.0, 1.0, 3.0};
  A->appendRow(3, index1, val1);
  A->appendRow(3, index2, val2);
}

void generate_conic_constraint(CoinPackedMatrix *& A, CoinPackedVector *& b,
                               CoinPackedVector *& d, double & h) {
  fill_A(A);
  b = new CoinPackedVector();
  int b_indices[] = {0,1};
  double b_values[] = {1.0, 2.0};
  b->setVector(2, b_indices, b_values);
  d = new CoinPackedVector();
  int d_indices[] = {0,1,3};
  double d_values[] = {-1.0, -1.0, 1.0};
  d->setVector(3, d_indices, d_values);
  h = 2.0;
}

int main () {
  // generate conic constraints
  CoinPackedMatrix * A = 0;
  CoinPackedVector * b = 0;
  CoinPackedVector * d = 0;
  double h = 0;
  generate_conic_constraint(A, b, d, h);
  A->dumpMatrix();
  // generate linear constraints
  CoinPackedMatrix * matrix = 0;
  double * col_lb = 0;
  double * col_ub = 0;
  double * obj = 0;
  double * row_lb = 0;
  double * row_ub = 0;
  generate_linear_constraints(matrix, col_lb, col_ub, obj, row_lb, row_ub);
  // create solver and load linear constraints
  OsiConicSolverInterface * solver = new ColaModel();
  // just add cols to the model
  solver->addCol(0, 0, 0, 0.0, solver->getInfinity(), 1.0, std::string("x1"));
  solver->addCol(0, 0, 0, 0.0, solver->getInfinity(), 0.0, std::string("x2"));
  solver->addCol(0, 0, 0, 0.0, solver->getInfinity(), 1.0, std::string("x3"));
  solver->addCol(0, 0, 0, 0.0, solver->getInfinity(), 0.0, std::string("x4"));
  //solver->loadProblem(*matrix, col_lb, NULL, obj, row_lb, row_ub);
  // add conic constraints
  solver->addConicConstraint(A, b, d, h);
  // todo(aykut) implement writing problem to mps files
  //solver->writeMps("ex");
  // solve problem
  solver->initialSolve();
  // print status
  std::cout << "is optimal: " << solver->isProvenOptimal() << std::endl;
  // print solution
  ColaModel * cm = dynamic_cast<ColaModel*>(solver);
  cm->print_stats();
  double const * sol = solver->getColSolution();
  for (int i=0; i<4; ++i) {
    std::cout << "x" << i << " " << sol[i] << std::endl;
  }
  delete A;
  delete b;
  delete d;
  delete matrix;
  delete solver;
  return 0;
}
