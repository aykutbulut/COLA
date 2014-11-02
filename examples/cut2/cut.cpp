#include <ColaModel.hpp>
// OsiConic headers
#include <OsiConicSolverInterface.hpp>
#include <OsiConicCuts.hpp>
// STDLIB headers
#include <iomanip>
#include <numeric>

#define WIDTH 10
#define SPACE "          "
#define LINE  "__________"

void print1(OsiConicSolverInterface * si);
void print_optimality_info(OsiConicSolverInterface * si);
void print_solution(OsiConicSolverInterface * si);
void print_basis_info(OsiConicSolverInterface * si);
void print_reduced_cost(OsiConicSolverInterface * si);
void print_tableau_info(OsiConicSolverInterface * si);
void print_Binv(OsiConicSolverInterface * si);
// cut method
OsiConicCuts * MIR_cut(OsiConicSolverInterface * si);

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
  // print tableau
  print1(si);
  // time to cut
  OsiConicCuts * cuts = MIR_cut(si);
  std::cout << "After cut" << std::endl;
  //print1(si);
  //delete cuts;
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
  std::cout << SPACE;
  for (int i=0; i<num_cols; ++i)
    std::cout << std::setw(WIDTH) << std::setprecision(3) << sol[i];
  // print slack information
  // compute row and then the difference with the left hand side
  double const * act = si->getRowActivity();
  char const * sense = si->getRowSense();
  double const * rhs = si->getRightHandSide();
  for (int i=0; i<num_rows; ++i)
    std::cout << std::setw(WIDTH) << rhs[i]-act[i];
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
  std::cout << SPACE;
  for (int i=0; i<num_cols; ++i)
    std::cout << std::setw(WIDTH) << cstat[i];
  std::cout << std::endl;
  delete[] rstat;
  delete[] cstat;
}

void print_reduced_cost(OsiConicSolverInterface * si) {
  int num_cols = si->getNumCols();
  const double * rc = si->getReducedCost();
  std::cout << "RCost  " << "|";
  std::cout << SPACE;
  for (int i=0; i<num_cols; ++i)
    std::cout << std::setw(WIDTH) << rc[i];
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
  double const * row_lb = si->getRowLower();
  double const * row_ub = si->getRowUpper();
  for (int i=0; i<num_rows; ++i) {
    si->getBInvARow(i, row, slack);
    double * B_inv_row = new double[num_rows]();
    si->getBInvRow(i, B_inv_row);
    std::cout << "Row " << std::setw(3) << i << "|";
    // print Binv_i times row lower bound
    double lb_i = std::inner_product(row_lb, row_lb+num_rows, B_inv_row, 0.0);
    std::cout << std::setw(WIDTH) << std::setprecision(3) << lb_i;
    // print BinvA
    for (int j=0; j<num_cols; ++j) {
      std::cout << std::setw(WIDTH) << row[j];
    }
    for (int j=0; j<num_rows; ++j) {
      std::cout << std::setw(WIDTH) << slack[j];
    }
    // print Binv_i times row upper bound
    double ub_i = std::inner_product(row_ub, row_ub+num_rows, B_inv_row, 0.0);
    std::cout << std::setw(WIDTH) << std::setprecision(3) << ub_i;
    // print Binv b
    double val = std::inner_product(rhs, rhs+num_rows, B_inv_row, 0.0);
    delete[] B_inv_row;
    std::cout << "|" << std::setw(14) << val;
    // print row activity
    std::cout << "|" << std::setw(WIDTH) << act[i];
    // pprint row sense
    std::cout << "|" << std::setw(1) << sense[i];
    // print rhs
    std::cout << "|" << std::setw(WIDTH) << rhs[i];
    // print row status information
    std::cout << "|" << std::setw(1) << rstat[i];
    std::cout << std::endl;
  }
  delete[] cstat;
  delete[] rstat;
  delete[] row;
  delete[] slack;
}


void print_Binv(OsiConicSolverInterface * si) {
  int num_rows = si->getNumRows();
  // Binv columns
  double ** cols;
  cols = new double*[num_rows];
  for (int i=0; i<num_rows; ++i) {
    cols[i] = new double[num_rows];
    si->getBInvCol(i, cols[i]);
  }
  int * basis;
  basis = new int[num_rows];
  si->getBasics(basis);
  // first print variables
  std::cout << "Binv" << std::endl;
  for (int i=0; i<num_rows; ++i) {
    std::cout << std::setw(WIDTH) << basis[i];
  }
  std::cout << std::endl;
  for (int i=0; i<num_rows; ++i) {
    std::cout << LINE;
  }
  std::cout << std::endl;
  // print Binv
  for (int i=0; i<num_rows; ++i) {
    // print ith row
    for (int j=0; j<num_rows; ++j) {
      std::cout << std::setw(WIDTH) << cols[j][i];
    }
    std::cout << std::endl;
  }
}


void print1(OsiConicSolverInterface * si) {
  print_optimality_info(si);
  print_solution(si);
  print_basis_info(si);
  print_reduced_cost(si);
  print_tableau_info(si);
  print_Binv(si);
}


// MIT cut method, creates MIR cut for a given cone.
// 1. first get simplex tableau
// 2. Locate which variable in cone is a basic variable.
// 3. Write basic variable in terms of the others.
//
//
//
//
//
// Keep a list of cuts that MIR is generted. We may not need to
// insert t again.
//

OsiConicCuts * MIR_cut(OsiConicSolverInterface * si) {
  OsiConicCuts * cuts = new OsiConicCuts();
  // Get current basis status, if cstat==1 then basic
  int num_rows = si->getNumRows();
  int num_cols = si->getNumCols();
  int * basis_index = new int[num_rows];
  si->getBasics(basis_index);
  // decide on cut generating cone, for now we use first
  int cut_cone = 0;
  OsiConeType type;
  int num_members;
  int * members;
  si->getConicConstraint(cut_cone, type, num_members, members);
  if (type!=OSI_QUAD) {
    std::cerr << "We support cones in canonical form only." << std::endl;
    throw "";
  }
  // decide basic variable we will replace, choose first available
  int ind = -1;
  std::set<int> basis_set(basis_index, basis_index+num_rows);
  for (int i=1; i<num_members; ++i) {
    if (basis_set.find(members[i])!=basis_set.end()) {
      ind = i;
      break;
    }
  }
  if (ind==-1) {
    std::cerr << "A suitable member have not been found." << std::endl;
    throw "";
  }
  // add rows that corresponds to polyhedral second order cone
  // add slack variables s_N. For this we add the following rows
  // add rows corresponding to nonbasic slack variables
  // == get nonbasic rows
  CoinPackedMatrix const * mat = si->getMatrixByRow();
  for (int i=0; i<num_rows; ++i) {
    if (basis_set.find(i+num_cols)!=basis_set.end()) {
      // slack i is nonbasic. add row i to the problem again.
      double lb = si->getRowLower()[i];
      double ub = si->getRowUpper()[i];
      int num_elem  = mat->getVectorLast(i)-mat->getVectorFirst(i);
      //int num_elem = mat->getVectorLengths()[i];
      int * cols = new int[num_elem];
      double * value = new double[num_elem];
      int k=0;
      for (int j=mat->getVectorFirst(i); j<mat->getVectorLast(i); ++j) {
	cols[k] = mat->getIndices()[j];
	value[k] = mat->getElements()[j];
	k++;
      }
      si->addRow(num_elem, cols, value, lb, ub);
      delete[] cols;
      delete[] value;
    }
  }
  // add slacks s_N

  // add t_i rows

  delete[] basis_index;
  return cuts;
}
