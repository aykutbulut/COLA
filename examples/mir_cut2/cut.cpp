#include <ColaModel.hpp>
// OsiConic headers
#include <OsiConicSolverInterface.hpp>
#include <OsiConicCuts.hpp>
// STDLIB headers
#include <iomanip>
#include <numeric>
#include <cmath>

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
double phi(double a, double alpha);

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
  //std::cout << "After cut" << std::endl;
  //print1(si);
  dynamic_cast<ColaModel*>(si)->report_feasibility();
  double const * sol = si->getColSolution();
  for (int i=0; i<si->getNumCols(); ++i) {
    std::cout << std::setw(WIDTH+5) << sol[i];
  }
  std::cout << std::endl;
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

// ==> when we say variables it only corresponds to explicit variables in the model
// excluding slacks inserted by clp. When we refer slack variables introduced
// by clp we use term slack variables.
// ==> Slack variables can also be enumerated by enumerating the rows, a unique
// slack variable for each row. In for loops I use this to enumerate slacks
// sometimes, be aware. row i corresponds to slack variable index by i+num_cols
// in the tableau.
// ==> Clp inserts a slack for each row.
// ==> in tableau indices from 0 to num_cols-1 correspond to variables and from
// num_cols to num_cols+num_rows-1 correspond to slack variables.

OsiConicCuts * MIR_cut(OsiConicSolverInterface * si) {
  si->writeMps("mir_before_cut");
  OsiConicCuts * cuts = new OsiConicCuts();
  // Get current basis status, if cstat==1 then basic
  int num_rows = si->getNumRows();
  int num_cols = si->getNumCols();
  int num_basic = num_rows;
  int num_nonbasic = num_cols;
  int * B_index = new int[num_basic];
  si->getBasics(B_index);
  // set of basic variables
  std::set<int> B_set(B_index, B_index+num_rows);
  int num_xB=0;
  for (int i=0; i<num_basic; ++i) {
    if (B_index[i]<num_cols) {
      num_xB++;
    }
  }
  int * xB = new int[num_xB];
  int num_sB = num_basic-num_xB;
  int * sB = new int[num_sB];
  for (int i=0,k1=0,k2=0; i<num_basic; ++i) {
    if (B_index[i]<num_cols) {
      xB[k1] = B_index[i];
      k1++;
    }
    else {
      sB[k2] = B_index[i];
      k2++;
    }
  }
  std::set<int> xB_set(xB, xB+num_xB);
  std::set<int> sB_set(sB, sB+num_sB);
  int num_xN = num_cols-num_xB;
  int * xN = new int[num_xN];
  for (int i=0,k=0; i<num_cols; ++i) {
    if (xB_set.find(i)==xB_set.end()) {
      xN[k]=i;
      k++;
    }
  }
  int num_sN = num_nonbasic-num_xN;
  int * sN = new int[num_sN];
  for (int i=0,k=0; i<num_rows; ++i) {
    if (sB_set.find(i+num_cols)==sB_set.end()) {
      sN[k] = i+num_cols;
      k++;
    }
  }
  // ==== decide on cut generating cone, for now we use first
  int cut_cone = 0;
  OsiConeType cut_cone_type;
  int cut_cone_size;
  int * cut_cone_members;
  si->getConicConstraint(cut_cone, cut_cone_type, cut_cone_size,
			 cut_cone_members);
  if (cut_cone_type!=OSI_QUAD) {
    std::cerr << "We support cones in canonical form only." << std::endl;
    throw "";
  }
  // ==== end of decide on cut generating cone, for now we use first
  // ==== decide basic variable we will replace, choose first available
  // for now we will use only one variable.
  int ind = -1;
  for (int i=1; i<cut_cone_size; ++i) {
    if (B_set.find(cut_cone_members[i])!=B_set.end()) {
      ind = cut_cone_members[i];
      break;
    }
  }
  if (ind==-1) {
    std::cerr << "A suitable member have not been found." << std::endl;
    throw "";
  }
  // check if any of the nonbasics is integer
  int flag=0;
  for (int i=0; i<num_xN; ++i) {
    if (si->isInteger(xN[i]))
      flag=1;
  }
  if (flag==0) {
    std::cerr << "There is no integer nonbasic, we can not generate cut."
	      << std::endl;
    throw "";
  }
  // determine row corresponding to basic variable ind
  for (int i=0; i<num_basic; ++i) {
    if (ind==B_index[i]) {
      ind=i;
      break;
    }
  }
  // ==== decide basic variable we will replace, choose first available
  // retreive ind related tableau information to be used.
  double * BInvA_x_i = new double[num_cols];
  double * BInvA_s_i = new double[num_rows];
  double * BInvN_x_i = new double[num_xN];
  double * BInvN_s_i = new double[num_sN];
  si->getBInvARow(ind, BInvA_x_i, BInvA_s_i);
  // get the columns corresponding to nonbasics.
  for (int i=0; i<num_xN; ++i) {
    BInvN_x_i[i] = BInvA_x_i[xN[i]];
  }
  for (int i=0; i<num_sN; ++i) {
    BInvN_s_i[i] = BInvA_s_i[sN[i]-num_cols];
  }
  // print
  std::cout << "BInvN_x_i" << std::endl;
  for (int i=0; i<num_xN; ++i) {
    std::cout << std::setw(WIDTH) << BInvN_x_i[i];
  }
  std::cout << std::endl;
  std::cout << "BInvN_s_i" << std::endl;
  for (int i=0; i<num_sN; ++i) {
    std::cout << std::setw(WIDTH) << BInvN_s_i[i];
  }
  std::cout << std::endl;
  // end of print
  double * BInv_i = new double[num_basic];
  si->getBInvRow(ind, BInv_i);
  double const * rhs = si->getRightHandSide();
  std::cout << "BInv_i" << std::endl;
  for (int i=0; i<num_basic; ++i) {
    std::cout << std::setw(WIDTH) << BInv_i[i];
  }
  std::cout << std::endl;
  // there should be at least one integer nonbasic with nonzero coef.
  flag=0;
  for (int i=0; i<num_xN; ++i) {
    if (si->isInteger(xN[i]) and BInvN_x_i[i]) {
      if (xN[i]==cut_cone_members[0])
	continue;
      flag=1;
    }
  }
  if (flag==0) {
    std::cerr << "There should be at least one integer nonbasic with nonzero coef!"
	      << std::endl;
    throw "";
  }
  double BInvRhs_i = std::inner_product(rhs, rhs+num_rows, BInv_i, 0.0);
  std::cout << "BInvRhs_i: " << std::endl;
  std::cout << BInvRhs_i << std::endl;
  // add slack variables s_N. we do this since we can not insert a row to the
  // simplex tableau directly. We add new cols that correspond to s_N
  // and rows that defines s_N
  int * sN_index = new int[num_sN];
  CoinPackedMatrix const * mat;
  for (int i=0, k=0; i<num_rows; ++i) {
    if (sB_set.find(i+num_cols)==sB_set.end()) {
      std::cout << "Adding row for nonbasic slack " << i << std::endl;
      // slack i is nonbasic. add row i to the problem again.
      // we only need to add nonzero coefficient ones but we
      // add all for now.
      mat = si->getMatrixByRow();
      double const * rhs_local = si->getRightHandSide();
      //mat->dumpMatrix();
      int num_elem  = mat->getVectorLast(i)-mat->getVectorFirst(i);
      //int num_elem = mat->getVectorLengths()[i];
      int * cols = new int[num_elem];
      double * value = new double[num_elem];
      int first = mat->getVectorFirst(i);
      int last = mat->getVectorLast(i);
      for (int j=first; j<last; ++j) {
	cols[j-first] = mat->getIndices()[j];
	value[j-first] = mat->getElements()[j];
      }
      si->addRow(num_elem, cols, value, rhs[i], rhs[i]);
      delete[] cols;
      delete[] value;
      // now add column s_N_i to this row
      int rows[] = {si->getNumRows()-1};
      double row_value[] = {1.0};
      si->addCol(1, rows, row_value, -si->getInfinity(),
		 si->getInfinity(), 0.0);
      sN_index[k] = k+num_cols;
      k++;
    }
  }
  // now we have all nonbasic slack as variables and we can use them, they
  // range from num_cols to si->getNumCols()
  // time to replace cone member ind with its nonbasic variable representation
  // in polyhedral second order conic representation, i is ind
  // t_i > |(x_B)i|
  // t_i > |(B^-1b)_i - (B^-1N)_i (x_N s_N)|
  // we add following two rows
  // t_i + (B^-1N)_i (x_N s_N) > (B^-1b)_i, (B^-1N)_i
  // t_i - (B^-1N)_i (x_N s_N) > -(B^-1b)_i
  int * e1 = new int[num_nonbasic]();
  double * v1 = new double[num_nonbasic]();
  std::copy(xN, xN+num_xN, e1);
  std::copy(sN_index, sN_index+num_sN, e1+num_xN);
  std::copy(BInvN_x_i, BInvN_x_i+num_xN, v1);
  std::copy(BInvN_s_i, BInvN_s_i+num_sN, v1+num_xN);
  si->addRow(num_nonbasic, e1, v1, BInvRhs_i, si->getInfinity());
  double * neg_v1 = new double[num_nonbasic]();
  for (int i=0; i<num_nonbasic; ++i) {
    neg_v1[i] = -v1[i];
  }
  si->addRow(num_nonbasic, e1, neg_v1, -BInvRhs_i, si->getInfinity());
  delete[] e1;
  delete[] v1;
  delete[] neg_v1;
  // add column t_i, it is in last two rows
  int e2[] = {si->getNumRows()-2, si->getNumRows()-1};
  double v2[] = {1.0, 1.0};
  si->addCol(2, e2, v2, 0.0, si->getInfinity(), 0.0);
  // add cone with t_i as a member
  // now time to add actual MIR cut
  double alpha = 1.5;
  double b_over_alpha = -BInvRhs_i/alpha;
  int b_over_alpha_lower = floor(b_over_alpha);
  double f_alpha = b_over_alpha - b_over_alpha_lower;
  //double b_over_alpha = 0.0;
  //double * a_over_alpha = 0;
  // find integer x_N.
  int * xN_int = 0;
  int num_xN_int = 0;
  for (int i=0; i<num_cols; ++i) {
    if (xB_set.find(i)==xB_set.end()) {
      if (si->isInteger(i))
	num_xN_int++;
    }
  }
  xN_int = new int[num_xN_int];
  for (int i=0, k=0; i<num_cols; ++i) {
    if (xB_set.find(i)==xB_set.end()) {
      if (si->isInteger(i)) {
	xN_int[k] = i;
	k++;
      }
    }
  }
  std::set<int> xN_int_set(xN_int, xN_int+num_xN_int);
  int * e3 = new int[num_nonbasic+1];
  std::copy(xN, xN+num_xN, e3);
  std::copy(sN_index, sN_index+num_sN, e3+num_xN);
  e3[num_nonbasic] = num_cols+num_sN;
  double * v3 = new double[num_nonbasic+1];
  for (int i=0; i<num_xN; ++i) {
    if (xN_int_set.find(xN[i])!=xN_int_set.end()) {
      v3[i] = phi(-BInvN_x_i[i]/alpha, f_alpha);
    }
    else {
      v3[i] = -1.0/abs(alpha);
    }
  }
  for (int i=num_xN; i<num_xN+num_sN; ++i) {
    v3[i] = -1.0/abs(alpha);
  }
  v3[num_nonbasic] = -1.0/abs(alpha);
  double phi_b_over_alpha = phi(b_over_alpha, f_alpha);
  // add mir cut itself
  si->addRow(num_nonbasic+1, e3, v3, -si->getInfinity(), phi_b_over_alpha);
  delete[] B_index;
  delete[] xB;
  delete[] xN;
  delete[] sB;
  delete[] sN;
  delete[] BInvN_x_i;
  delete[] BInvN_s_i;
  delete[] BInv_i;
  delete[] sN_index;
  // modify cone of the model, replace column 1 (x_2) with column 7 (t_1)
  // since we added three nonbasic slacks (6,14,15) and then t_1
  si->removeConicConstraint(0);
  int * mem = new int[3];
  mem[0]=0;
  mem[1]=si->getNumCols()-1;
  mem[2]=2;
  si->addConicConstraint(OSI_QUAD, 3, mem);
  // resolve
  double old_obj = si->getObjValue();
  std::cout << "Before cut: " << std::setprecision(10) << old_obj <<std::endl;
  // print cut
  std::cout << "Cut: " << std::endl;
  std::cout << "Members: ";
  for (int i=0; i<num_nonbasic+1; ++i) {
    std::cout << std::setw(WIDTH) << e3[i];
  }
  std::cout << std::endl;
  std::cout << "Values:  ";
  for (int i=0; i<num_nonbasic+1; ++i) {
    std::cout << std::setw(WIDTH) << v3[i];
  }
  std::cout << "     >= " << phi_b_over_alpha << std::endl;
  delete[] e3;
  delete[] v3;

  // write mps file
  si->writeMps("mir_after_cut");
  si->resolve();
  std::cout << "After cut: " << si->getObjValue() << std::endl;
  return cuts;
}

double phi(double a, double f) {
  int n= floor(a);
  double value;
  if (a<n+f) {
    value = (1-2*f)*n-(a-n);
  }
  else {
    value = (1-2*f)*n+(a-n)-2*f;
  }
  return value;
}


