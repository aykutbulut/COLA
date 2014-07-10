#include "ColaModel.hpp"
#include "ConicConstraints.hpp"
#include "Separate.hpp"
#include "CoinMpsIO.hpp"
#include "Cut.hpp"

//#include <OsiMskSolverInterface.hpp>
#include <OsiClpSolverInterface.hpp>
#include <CoinPackedMatrix.hpp>

#include <cstring>
#include <iomanip>
#include <numeric>
#include <cmath>

#define TINY_EPS 1e-15

ColaModel::ColaModel() : cc_(NULL) {
  options_ = new Options();
  total_num_supports_ = 0;
  total_num_cuts_ = 0;
  num_supports_ = 0;
  num_cuts_ = 0;
  // this is the default beavior, user can change this using options
  setHintParam(OsiDoReducePrint,true,OsiHintTry);
  cc_ = new ConicConstraints();
}

ColaModel::ColaModel(char * data_file) : cc_(NULL) {
  options_ = new Options();
  total_num_supports_ = 0;
  total_num_cuts_ = 0;
  num_supports_ = 0;
  num_cuts_ = 0;
  // this is the default beavior, user can change this using options
  setHintParam(OsiDoReducePrint,true,OsiHintTry);
  cc_ = new ConicConstraints();
  read(data_file);
}

// overrite clone function of OsiClpSolverInterface
ColaModel * ColaModel::clone (bool copyData) const {
  ColaModel * new_model = new ColaModel();
  // copy solver
  // copy conic constraints
  new_model->setCC(cc_->clone());
  // copy options
  new_model->setOptions(options_->clone());
  // set cut and support statistics
  new_model->set_num_cuts(num_cuts_);
  new_model->set_num_supports(num_supports_);
  new_model->set_total_num_cuts(total_num_cuts_);
  new_model->set_total_num_supports(total_num_supports_);
  return new_model;
}

ColaModel::~ColaModel() {
  delete cc_;
  if(num_cuts_) {
    delete[] num_cuts_;
  }
  if (num_supports_) {
    delete[] num_supports_;
  }
}

// this read uses coinutils, ie coiniomps. Thanks to it  we do not need mosek
// to read conic problems, and now we can use CLP to solve our LP problems.
// this function is based on readconic.cpp example in Osi.
void ColaModel::read(const char * data_file) {
  CoinMpsIO m_MpsData;
  int nOfSOS;
  CoinSet ** SOS = NULL;
  int status = m_MpsData.readMps(data_file, "", nOfSOS, SOS );
  if (nOfSOS) {
    throw "Input file has SOS section!";
  }
  delete [] SOS;
  assert (!status);
  int nOfCones;
  int * coneStart = NULL;
  int * coneIdx = NULL;
  int * coneType = NULL;
  status = m_MpsData.readConicMps(NULL, coneStart, coneIdx, coneType, nOfCones);
  assert (!status);
  int * members;
  for (int i=0; i<nOfCones; ++i) {
    if (coneType[i]!=1 and coneType[i]!=2) {
      throw "Invalid cone type!";
    }
    int num_members = coneStart[i+1]-coneStart[i];
    if (coneType[i]==2 and num_members<3) {
      throw "Rotated cones should have at least 3 members!";
    }
    // get members
    members = new int[num_members];
    int k=0;
    for (int j=coneStart[i]; j<coneStart[i+1]; ++j) {
      members[k] = coneIdx[j];
      k++;
    }
    ConeType type;
    if (coneType[i]==1) {
      type = QUAD;
    }
    else if (coneType[i]==2) {
      type = RQUAD;
    }
    cc_-> add_cone(num_members, members, type);
    delete[] members;
  }
  // check log level and print ccordingly
  if (nOfCones) {
    printf("Conic section has %d cones\n",nOfCones);
    for (int iCone=0;iCone<nOfCones;iCone++) {
      printf("Cone %d has %d entries (type %d) ",iCone,coneStart[iCone+1]-coneStart[iCone],
	     coneType[iCone]);
      for (int j=coneStart[iCone];j<coneStart[iCone+1];j++)
	printf("%d ",coneIdx[j]);
      printf("\n");
    }
  }
  delete [] coneStart;
  delete [] coneIdx;
  delete [] coneType;
  // load problem
  const CoinPackedMatrix * matrix = m_MpsData.getMatrixByCol();
  const double * collb = m_MpsData.getColLower();
  const double * colub = m_MpsData.getColUpper();
  const double * obj = m_MpsData.getObjCoefficients();
  const double * rowlb = m_MpsData.getRowLower();
  const double * rowub = m_MpsData.getRowUpper();
  loadProblem(*matrix, collb, colub, obj, rowlb, rowub);
  // set row and column names
  // todo(aykut) names are assumed to be less than 255 characters
  int name_len = 255;
  int numcols=m_MpsData.getNumCols();
  int numrows=m_MpsData.getNumRows();
  for (int i=0; i<numrows; ++i) {
    setRowName(i, m_MpsData.rowName(i));
  }
  for (int i=0; i<numcols; ++i) {
    setColName(i, m_MpsData.columnName(i));
  }
  // set variable types
  for (int i=0; i<numcols; ++i) {
    if (m_MpsData.isInteger(i)) {
      setInteger(i);
    }
  }
  // allocate memory for cut and support statistics
  num_cuts_ = new int[nOfCones]();
  num_supports_ = new int[nOfCones]();
}

void ColaModel::print_stats() const {
  // how many separ problem solved
  // separation for each cone
  int num_cones = cc_->num_cones();
  std::cout << "Total Unboundedness supports: " << total_num_supports_ << std::endl;
  for (int i=0; i<num_cones; ++i) {
    std::cout << "     Supports for cone " << std::setw(5) << i
	      << std::setw(6) << num_supports_[i] << std::endl;
  }
  std::cout << "Total Seperating hyperplanes: " << total_num_cuts_ << std::endl;
  for (int i=0; i<num_cones; ++i) {
    std::cout << "     Seperating hyperplane for cone " << std::setw(5) << i
	      << std::setw(6) << num_cuts_[i] << std::endl;
  }
}

void ColaModel::print_solution() const {
  std::cout << "Solution is" << std::endl;
  const double * sol = getColSolution();
  for (int i=0; i < getNumCols(); ++i) {
    std::cout << std::setw(5) << i << std::setw(15) << sol[i] << std::endl;
  }
  std::cout << "Objective Value " << getObjValue() << std::endl;
}

ProblemStatus ColaModel::solve() {
  // solve linear problem
  bool feasible = false;
  Separate * sep;
  // add nonnegativity of leading variables
  int num_cones = cc_->num_cones();
  for (int i=0; i<num_cones; ++i) {
    if(cc_->type(i)==QUAD) {
      setColLower(cc_->cone_members(i)[0], 0.0);
    }
    else {
      setColLower(cc_->cone_members(i)[0], 0.0);
      setColLower(cc_->cone_members(i)[1], 0.0);
    }
  }
  // ===== End Of adding nonnegativity of leading variables
  //writeMps("initial", "mps");
  OsiClpSolverInterface::initialSolve();
  // check problem status
  problem_status();
  if ((soco_status_!=OPTIMAL) && (soco_status_!=DUAL_INFEASIBLE))
    return soco_status_;
  // if problem is unbounded add supporting hyperplanes for each cone using
  // objective coefficients.
  // todo(aykut) I guess SOCO is unbounded if LP is unbounded after all these
  // supporting hyperplanes. Check if this is true theoretically.
  if (soco_status_==DUAL_INFEASIBLE) {
    std::cout << "Cola: Problem without conic constraints is unbounded. Adding"
      " supporting hyperplanes to resolve..."
	      << std::endl;
    while (soco_status_==DUAL_INFEASIBLE) {
      // get one ray
      // todo(aykut) for now we get only one ray
      std::vector<double*> rays = getPrimalRays(1);
      const double * vec = 0;
      if (!rays.empty() and rays[0]!=0) {
	vec = rays[0];
      }
      else {
	// check if primal is infeasible, then the problem is infeasible
	if (isProvenPrimalInfeasible()) {
	  std::cout << "Cola: Both LP primal and dual is infeasible."
		    << "Cola: This means cone problem is infeasible."
		    << std::endl
		    << "Cola: Terminating...";
	}
	else {
	  std::cout << "Cola: Warning! "
		    << "LP is dual infeasible but solver did not return a "
	    "direction of unboundedness." << std::endl
		    << "Cola: Trying to generate supports using objective "
	    "function coefficients..." << std::endl;
	  vec = getObjCoefficients();
	}
      }
      sep = new Separate(cc_, vec, getNumCols(), options_);
      // returns true if direction is feasible for all cones.
      feasible = sep->is_feasible();
      if (feasible) {
	// primal ray is feasible for all cone constraints,
	// problem is unbounded
	return DUAL_INFEASIBLE;
      }
      else {
	// Add all the cuts generated
	std::vector<Cut*>::const_iterator it;
	for(it=sep->cuts()->begin(); it!=sep->cuts()->end(); ++it) {
	  addRow((**it).size(), (**it).index(), (**it).coef(),
			  -getInfinity(), (**it).rhs());
	  total_num_supports_++;
	  num_supports_[(*it)->cut_generating_cone()]++;
	}
      }
      // todo(aykut) delete all rays not just first one.
      if (!rays.empty()) {
	delete[] rays[0];
	rays.empty();
      }
      delete sep;
      OsiClpSolverInterface::resolve();
      // update soco_status_
      problem_status();
    }
  }
  // it means it is not optimal after resolving unbounded directions
  if (soco_status_!=OPTIMAL) {
    std::cout << "Cola: Problem status is not optimal after adding "
      "supporting hyperplanes."
	      << std::endl;
    std::cout << "Cola: Terminating..." << std::endl;
    return soco_status_;
  }
  // when we reached here, it means the problem is not unbounded anymore.
  // todo(aykut) pick the cone and all those logs and stuff. it can all go to separate.
  sep = new Separate(cc_, getColSolution(), getNumCols(), options_);
  // returns true if given point is feasible for all cones.
  feasible = sep->is_feasible();
  while(!feasible) {
    // add all the cuts generated
    std::vector<Cut*>::const_iterator it;
    for(it=sep->cuts()->begin(); it!=sep->cuts()->end(); ++it) {
      addRow((**it).size(), (**it).index(), (**it).coef(),
		      -getInfinity(), (**it).rhs());
      total_num_cuts_++;
      num_cuts_[(**it).cut_generating_cone()]++;
    }
    // resolve the problem
    OsiClpSolverInterface::resolve();
    // update problem status
    problem_status();
    // check problem status
    if (soco_status_!=OPTIMAL)
      break;
    delete sep;
    // chec if the new solution is feasible for cones
    sep = new Separate(cc_, getColSolution(), getNumCols(), options_);
    feasible = sep->is_feasible();
  }
  // update problem status
  problem_status();
  return soco_status_;
}

ConicConstraints * ColaModel::get_conic_constraints() {
  return cc_;
}

Options * ColaModel::options() {
  return options_;
}


void ColaModel::setCC(ConicConstraints * cc) {
  cc_ = cc;
}

void ColaModel::setOptions(Options * opt) {
  options_ = opt;
}

void ColaModel::report_feasibility() const {
  int num_cones = cc_->num_cones();
  double * par_sol;
  double lhs = 0.0;
  double lhs_real = 0.0;
  std::cout << "Conic Constraints feasibility report" << std::endl;
  std::cout << std::setw(5) << std::left << "Cone"
            << std::setw(15) << std::left << "lhs"
            << std::setw(15) << std::left << "lhs_real"
            << std::endl;
  for (int i=0; i<num_cones; ++i) {
    int cone_size = cc_->cone_size(i);
    const int * members = cc_->cone_members(i);
    par_sol = new double[cone_size];
    const double * full_sol = getColSolution();
    for (int j=0; j<cone_size; ++j) {
      par_sol[j] = full_sol[members[j]];
    }
    if (cc_->type(i)==QUAD) {
      lhs = par_sol[0]*par_sol[0]
	- std::inner_product(par_sol+1, par_sol+cone_size, par_sol+1, 0.0);
      lhs_real = par_sol[0]
	-sqrt(std::inner_product(par_sol+1, par_sol+cone_size, par_sol+1, 0.0));
    }
    else if (cc_->type(i)==RQUAD) {
      lhs = 2.0*par_sol[0]*par_sol[1]
	- std::inner_product(par_sol+2, par_sol+cone_size, par_sol+2, 0.0);
      lhs_real = lhs;
    }
    std::cout << std::setw(5) << std::left << i
              << std::setw(15) << std::left << lhs
              << std::setw(15) << std::left << lhs_real
              << std::endl;
    delete[] par_sol;
  }
}

void ColaModel::initialSolve() {
  solve();
}

void ColaModel::resolve() {
  solve();
}

int ColaModel::readMps(const char * filename, const char * extension) {
  read(filename);
  return 0;
}

// returns problem status and updates status_
ProblemStatus ColaModel::problem_status() {
  if (isAbandoned()) {
    soco_status_ = ABANDONED;
  }
  else if (isProvenOptimal()) {
    soco_status_ = OPTIMAL;
  }
  else if (isProvenPrimalInfeasible()) {
    soco_status_ = PRIMAL_INFEASIBLE;
  }
  else if (isProvenDualInfeasible()) {
    soco_status_ = DUAL_INFEASIBLE;
  }
  else if (isPrimalObjectiveLimitReached()) {
    soco_status_= PRIMAL_OBJECTIVE_LIMIT_REACHED;
  }
  else if (isDualObjectiveLimitReached()) {
    soco_status_= DUAL_OBJECTIVE_LIMIT_REACHED;
  }
  else if (isIterationLimitReached()) {
    soco_status_= ITERATION_LIMIT_REACHED;
  }
  return soco_status_;
}

const int * ColaModel::num_cuts() const {
  return num_cuts_;
}

int ColaModel::num_cuts(int i) const {
  return num_cuts_[i];
}

const int * ColaModel::num_supports() const {
  return num_supports_;
}

int ColaModel::num_supports(int i) const {
  return num_supports_[i];
}

int ColaModel::total_num_cuts() const {
  return total_num_cuts_;
}

int ColaModel::total_num_supports() const {
  return total_num_supports_;
}

// set cut and support statistics
void ColaModel::set_num_cuts(const int * num_cuts) {
  int nc = cc_->num_cones();
  std::copy(num_cuts, num_cuts+nc, num_cuts_);
}

void ColaModel::set_num_supports(const int * num_supports) {
  int nc = cc_->num_cones();
  std::copy(num_supports, num_supports+nc, num_supports_);
}

void ColaModel::set_total_num_cuts(int tnc) {
  total_num_cuts_ = tnc;
}

void ColaModel::set_total_num_supports(int tns) {
  total_num_supports_ = tns;
}

