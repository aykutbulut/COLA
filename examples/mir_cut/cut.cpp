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

void choose_cut_var(OsiConicSolverInterface * si, int & cut_cone, int & cut_var, int & cut_row);
// cut method
//for now it adds cuts directly to solver
void add_MIR_cut(OsiConicSolverInterface * si, int cut_cone,
		 int cut_var, int cut_row);
double phi(double a, double alpha);

int main(int argc, char ** argv) {
  OsiConicSolverInterface * si = new ColaModel();
  // read
  si->readMps(argv[1]);
  // decide variable to generate cut, it should be a member of cone.
  // decide row to use, it should have a nonzero coef for cut generating
  // variable
  int cut_cone = -1;
  int cut_var = -1;
  int cut_row = -1;
  choose_cut_var(si, cut_cone, cut_var, cut_row);
  std::cout << "Cut generating cone: " << cut_cone << std::endl;
  std::cout << "Cut generating var : " << cut_var << std::endl;
  std::cout << "Cut generating row : " << cut_row << std::endl;
  // generte and add cuts
  add_MIR_cut(si, cut_cone, cut_var, cut_row);
  si->initialSolve();
  std::cout << "After cut objective: " << si->getObjValue() << std::endl;
  // print objective after cut
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

void add_MIR_cut(OsiConicSolverInterface * si, int cut_cone, int cut_var, int cut_row) {
  int num_rows = si->getNumRows();
  int num_cols = si->getNumCols();
  //si->writeMps("mir_before_cut");
  // add t rows
  CoinPackedMatrix const * mat;
  mat = si->getMatrixByRow();
  double const * rhs = si->getRightHandSide();
  int first = mat->getVectorFirst(cut_row);
  int last = mat->getVectorLast(cut_row);
  int num_elem  = last-first;
  int * cols = new int[num_elem];
  double * value = new double[num_elem];
  for (int j=first; j<last; ++j) {
    cols[j-first] = mat->getIndices()[j];
    value[j-first] = mat->getElements()[j];
    if (cols[j-first]==cut_var)
      value[j-first] = value[j-first]-1.0;
  }
  si->addRow(num_elem, cols, value, rhs[cut_row], si->getInfinity());
  double * neg_value = new double[num_elem];
  for (int i=0; i<num_elem; ++i) {
    neg_value[i] = -value[i];
  }
  si->addRow(num_elem, cols, neg_value, -rhs[cut_row], si->getInfinity());
  // add t columns
  int t_e[2] = {num_rows, num_rows+1};
  double t_v[2]  = {1.0, 1.0};
  si->addCol(2, t_e, t_v, 0.0, si->getInfinity(), 0.0);
  // add cut
  int * cut_e = new int[num_elem+1];
  double * cut_v = new double[num_elem+1];
  std::copy(cols, cols+num_elem, cut_e);
  // index of variable t
  cut_e[num_elem] = num_cols;
  double alpha = 1.5;
  double f_alpha = rhs[cut_row]/alpha - floor(rhs[cut_row]/alpha);
  for (int i=0; i<num_elem; ++i) {
    if (si->isInteger(cols[i])) {
      cut_v[i] = phi(value[i]/alpha, f_alpha);
    }
    else {
      cut_v[i] = -1.0/abs(alpha);
    }
  }
  cut_v[num_elem] = -1.0/abs(alpha);
  si->addRow(num_elem+1, cut_e, cut_v, -si->getInfinity(), phi(rhs[cut_row]/alpha, f_alpha));
  // modify cone
  si->writeMps("mir_after_cut");
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

void choose_cut_var(OsiConicSolverInterface * si, int & cut_cone,
		    int & cut_var, int & cut_row) {
  int num_cones = si->getNumCones();
  int num_cols = si->getNumCols();
  int num_rows = si->getNumRows();
  char const * row_sense = si->getRowSense();
  // pick the first cone to generate cuts
  cut_cone = 0;
  OsiLorentzConeType cut_cone_type;
  int cut_cone_size = -1;
  int * cut_cone_members = 0;
  si->getConicConstraint(cut_cone, cut_cone_type, cut_cone_size,
			 cut_cone_members);
  if (cut_cone_type!=OSI_QUAD) {
    std::cerr << "We support cones in canonical form only." << std::endl;
    throw "";
  }
  // pick first cont var in the cut cone
  for(int i=1; i<cut_cone_size; ++i) {
    if (!si->isInteger(cut_cone_members[i])) {
      cut_var = cut_cone_members[i];
      break;
    }
  }
  // pick first row that has nonzero coef for cut_var
  cut_row=-1;
  CoinPackedMatrix const * mat;
  mat = si->getMatrixByRow();
  //double const * rhs_local = si->getRightHandSide();
  for (int i=0; i<num_rows; ++i) {
    if (row_sense[i]!='E')
      continue;
    int first = mat->getVectorFirst(i);
    int last = mat->getVectorLast(i);
    int num_elem  = last-first;
    int * cols = new int[num_elem];
    double * value = new double[num_elem];
    // check if cut_var is in the row
    int flag = 0;
    for (int j=first; j<last; ++j) {
      cols[j-first] = mat->getIndices()[j];
      value[j-first] = mat->getElements()[j];
      if (cols[j-first]==cut_var and value[j-first]!=0.0) {
	flag=1;
	break;
      }
    }
    if (flag) {
      cut_row=i;
      break;
    }
  }
  if (cut_row==-1) {
    std::cerr << "We could not find a row that has variable " << cut_var << std::endl;
    throw "";
  }
  delete[] cut_cone_members;
}
