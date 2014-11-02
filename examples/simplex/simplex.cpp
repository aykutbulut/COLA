#include <ColaModel.hpp>
#include <OsiConicSolverInterface.hpp>
// STDLIB headers
#include <iomanip>
#include <numeric>

void print_optimality_info(OsiConicSolverInterface * si);
void print_solution(OsiConicSolverInterface * si);
void print_basis_info(OsiConicSolverInterface * si);
void print_reduced_cost(OsiConicSolverInterface * si);
void print_tableau_info(OsiConicSolverInterface * si);

int main(int argc, char ** argv) {
  OsiConicSolverInterface * si = new ColaModel();
  //OsiSolverInterface * si = new OsiClpSolverInterface();
  // querry whether simplex info available
  std::cout << "Can do simplex interface? " << si->canDoSimplexInterface()
	    << std::endl;
  std::cout << "Is basis is available? " << si->basisIsAvailable()
	    << std::endl;
  std::cout << "Enabling simplex interface... " << std::endl;
  si->enableSimplexInterface(1);
  // Tells solver that calls to getBInv etc are about to take place.
  std::cout << "Enabling factorization... " << std::endl;
  si->enableFactorization();
  // read
  si->readMps(argv[1]);
  // solve
  si->initialSolve();
  si->enableSimplexInterface(1);
  // Tells solver that calls to getBInv etc are about to take place.
  std::cout << "Enabling factorization... " << std::endl;
  si->enableFactorization();
  print_optimality_info(si);
  print_solution(si);
  print_basis_info(si);
  print_reduced_cost(si);
  print_tableau_info(si);
  delete si;
  return 0;
}


void print_optimality_info(OsiConicSolverInterface * si) {
  std::cout << "Is Optimality proven? " << si->isProvenOptimal()
	    << std::endl;
  std::cout << "Is primal infeasible? " << si->isProvenPrimalInfeasible()
	    << std::endl;
  std::cout << "Is dual infeasible?   " << si->isProvenPrimalInfeasible()
	    << std::endl;
  std::cout << "Optimal value is " << si->getObjValue() << std::endl;
}

void print_solution(OsiConicSolverInterface * si) {
  int num_cols = si->getNumCols();
  int num_rows = si->getNumRows();
  double const * sol = si->getColSolution();
  std::cout << "Sol.   " << "|";
  for (int i=0; i<num_cols; ++i)
    std::cout << std::setw(14) << sol[i];
  // print slack information
  // compute row and then the difference with the left hand side
  double const * act = si->getRowActivity();
  char const * sense = si->getRowSense();
  double const * rhs = si->getRightHandSide();
  for (int i=0; i<num_rows; ++i)
    std::cout << std::setw(14) << rhs[i]-act[i];
  std::cout << std::endl;
}

void print_basis_info(OsiConicSolverInterface * si) {
  // we print column status here but row status will be printed
  // when we print the tableau
  int num_cols = si->getNumCols();
  int num_rows = si->getNumRows();
  int * cstat = new int[num_cols];
  int * rstat = new int[num_rows];
  si->getBasisStatus(cstat, rstat);
  std::cout << "Basis  " << "|";
  for (int i=0; i<num_cols; ++i)
    std::cout << std::setw(14) << cstat[i];
  std::cout << std::endl;
  delete[] rstat;
  delete[] cstat;
}

void print_reduced_cost(OsiConicSolverInterface * si) {
  int num_cols = si->getNumCols();
  const double * rc = si->getReducedCost();
  std::cout << "RCost  " << "|";
  for (int i=0; i<num_cols; ++i)
    std::cout << std::setw(14) << rc[i];
  std::cout << std::endl;
}

void print_tableau_info(OsiConicSolverInterface * si) {
  int num_cols = si->getNumCols();
  int num_rows = si->getNumRows();
  int * cstat = new int[num_cols];
  int * rstat = new int[num_rows];
  si->getBasisStatus(cstat, rstat);
  double * row = new double[num_cols];
  double * slack = new double[num_rows];
  double const * act = si->getRowActivity();
  char const * sense = si->getRowSense();
  double const * rhs = si->getRightHandSide();
  for (int i=0; i<num_rows; ++i) {
    si->getBInvARow(i, row, slack);
    // print BinvA
    std::cout << "Row " << std::setw(3) << i << "|";
    for (int j=0; j<num_cols; ++j) {
      std::cout << std::setw(14) << row[j];
    }
    for (int j=0; j<num_rows; ++j) {
      std::cout << std::setw(14) << slack[j];
    }
    // print Binv b
    double * B_inv_row = new double[num_cols]();
    si->getBInvRow(i, B_inv_row);
    double val = std::inner_product(rhs, rhs+num_rows, B_inv_row, 0.0);
    delete[] B_inv_row;
    std::cout << "|" << std::setw(14) << val;
    // print row activity
    std::cout << "|" << std::setw(14) << act[i];
    // pprint row sense
    std::cout << "|" << std::setw(1) << sense[i];
    // print rhs
    std::cout << "|" << std::setw(14) << rhs[i];
    // print row status information
    std::cout << "|" << std::setw(1) << rstat[i];
    std::cout << std::endl;
  }
  delete[] cstat;
  delete[] rstat;
  delete[] row;
  delete[] slack;
}
