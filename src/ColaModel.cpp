#include "ColaModel.hpp"
#include "ConicConstraints.hpp"
#include "Separate.hpp"
#include "CoinMpsIO.hpp"
#include "Cut.hpp"

#include <OsiClpSolverInterface.hpp>
#include <CoinPackedMatrix.hpp>
#include <CoinPackedVector.hpp>

#include <cstring>
#include <iomanip>
#include <numeric>
#include <cmath>

ColaModel::ColaModel() : cc_(NULL), OsiClpSolverInterface() {
  options_ = new Options();
  total_num_supports_ = 0;
  total_num_cuts_ = 0;
  num_supports_ = 0;
  num_cuts_ = 0;
  num_lp_solved_ = 0;
  // this is the default beavior, user can change this using options
  setHintParam(OsiDoReducePrint,true,OsiHintTry);
  cc_ = new ConicConstraints();
  // for unboundedness directions set option
  OsiClpSolverInterface::getModelPtr()->setMoreSpecialOptions(0);
}

ColaModel::ColaModel(char * data_file) : cc_(NULL), OsiClpSolverInterface() {
  options_ = new Options();
  total_num_supports_ = 0;
  total_num_cuts_ = 0;
  num_supports_ = 0;
  num_cuts_ = 0;
  num_lp_solved_ = 0;
  // this is the default beavior, user can change this using options
  setHintParam(OsiDoReducePrint,true,OsiHintTry);
  cc_ = new ConicConstraints();
  OsiConicSolverInterface::readMps(data_file);
  OsiClpSolverInterface::getModelPtr()->setMoreSpecialOptions(0);
}

// copy constructor
ColaModel::ColaModel(const ColaModel & other): OsiClpSolverInterface(other) {
  // copy conic constraints
  cc_ = other.get_conic_constraints()->clone();
  // copy options
  options_ = other.options()->clone();
  // copy problem status
  soco_status_ = other.problem_status();
  // copy number of lp solved
  num_lp_solved_ = other.num_lp_solved();
  // copy number of cuts generated
  int numCones = other.getNumCones();
  num_cuts_ = new int[numCones]();
  std::copy(other.num_cuts(), other.num_cuts()+numCones, num_cuts_);
  total_num_cuts_ = other.total_num_cuts();
  // copy number of supports
  num_supports_ = new int[numCones]();
  std::copy(other.num_supports(), other.num_supports() + numCones,
	    num_supports_);
  total_num_supports_ = other.total_num_supports();
}

// copy assignment operator
ColaModel & ColaModel::operator=(const ColaModel & rhs) {
  // copy conic constraints
  if (cc_)
    delete cc_;
  cc_ = rhs.get_conic_constraints()->clone();
  // copy options
  if (options_)
    delete options_;
  options_ = rhs.options()->clone();
  // copy problem status
  soco_status_ = rhs.problem_status();
  // copy number of lp solved
  num_lp_solved_ = rhs.num_lp_solved();
  int numCones = rhs.getNumCones();
  // copy number of cuts generated
  if (num_cuts_)
    delete num_cuts_;
  num_cuts_ = new int[numCones]();
  std::copy(rhs.num_cuts(), rhs.num_cuts()+numCones, num_cuts_);
  total_num_cuts_ = rhs.total_num_cuts();
  // copy number of supports
  if (num_supports_)
    delete num_supports_;
  num_supports_ = new int[numCones]();
  std::copy(rhs.num_supports(), rhs.num_supports() + numCones,
	    num_supports_);
  total_num_supports_ = rhs.total_num_supports();
  return *this;
}

// overrite clone function of OsiConicSolverInterface
OsiConicSolverInterface * ColaModel::clone (bool copyData) const {
  ColaModel * new_solver;
  if (copyData) {
    new_solver = new ColaModel(*this);
  }
  else {
    new_solver = new ColaModel();
  }
  return new_solver;
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


void ColaModel::addConicConstraint(OsiConeType type,
				  int numMembers,
				   const int * members) {
  ConeType t;
  if (type==OSI_QUAD)
    t = QUAD;
  else
    t = RQUAD;
  cc_-> add_cone(numMembers, members, t);
}

// this read uses coinutils, ie coiniomps. Thanks to it  we do not need mosek
// to read conic problems, and now we can use CLP to solve our LP problems.
// this function is based on readconic.cpp example in Osi.
// void ColaModel::read(const char * data_file) {
//   readMps(data_file);
//   // allocate memory for cut and support statistics
//   int nOfCones = getNumCones();
//   num_cuts_ = new int[nOfCones]();
//   num_supports_ = new int[nOfCones]();
// }

int ColaModel::readMps(const char * filename, const char * extension) {
  OsiConicSolverInterface::readMps(filename, extension);
  int nOfCones = getNumCones();
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
  std::cout << "Number of LPs solved: " << num_lp_solved_ << std::endl;
}

void ColaModel::print_solution() const {
  std::cout << "Solution is" << std::endl;
  const double * sol = getColSolution();
  for (int i=0; i < getNumCols(); ++i) {
    std::cout << std::setw(5) << i << std::setw(15) << sol[i] << std::endl;
  }
  std::cout << "Objective Value " << getObjValue() << std::endl;
}

ProblemStatus ColaModel::solve(bool resolve) {
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
  num_lp_solved_++;
  if (resolve==false) {
    // if this function is called from OsiConicSolverInterface::initialSolve then
    OsiClpSolverInterface::initialSolve();
  }
  else {
    // if it is called from OsiConicSolverInterface::resolve
    OsiClpSolverInterface::resolve();
  }
  // check problem status
  update_problem_status();
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
      // check if primal is infeasible, then the problem is infeasible
      if (isProvenPrimalInfeasible()) {
	// Both LP primal and dual is infeasible, conic problem is infeasible
	std::cout << "Cola: Conic problem is infeasible."
		  << std::endl
		  << "Cola: Terminating...";
      }

      // get one ray
      // todo(aykut) for now we get only one ray
      std::vector<double*> rays = getPrimalRays(1);
      const double * vec = 0;
      if (!rays.empty() and rays[0]!=0) {
	vec = rays[0];
      }
      else {
	std::cout << "Cola: Warning! "
		  << "LP is unbounded but solver did not return a "
	  "direction of unboundedness." << std::endl
		  << "Cola: Trying to generate supports using objective "
	  "function coefficients..." << std::endl;
	vec = getObjCoefficients();
      }
      // PRINT UNBDDNESS DIRECTION
      std::cout << "Unboundedness direction is " << std::endl  << "[";
      for (int i=0; i<getNumCols(); ++i) {
	std::cout << std::setw(10) << vec[i] << "; ";
      }
      std::cout << "]" << std::endl;
      // END OF DIRECTION PRINT
      sep = new Separate(cc_, vec, getNumCols(), options_);
      // returns true if direction is feasible for all cones.
      feasible = sep->is_feasible();
      if (feasible) {
	// primal ray is feasible for all cone constraints,
	// problem is unbounded
	delete sep;
	return DUAL_INFEASIBLE;
      }
      else {
	// Add all the cuts generated
	std::vector<Cut*>::const_iterator it;
	for(it=sep->cuts()->begin(); it!=sep->cuts()->end(); ++it) {
	  addRow((**it).size(), (**it).index(), (**it).coef(),
			  -getInfinity(), (**it).rhs());
	  // PRINT CUT
	  std::cout << "Cut " << std::endl  << "[";
	  for (int i=0; i<getNumCols(); ++i) {
	    std::cout << std::setw(10) << (**it).coef()[i] << "; ";
	  }
	  std::cout << "]" << std::endl;
	  // END OF CUT PRINT
	  total_num_supports_++;
	  num_supports_[(*it)->cut_generating_cone()]++;
	}
      }
      // todo(aykut) delete all rays not just first one.
      if (!rays.empty()) {
	for(int i=0; i<rays.size(); ++i)
	  delete[] rays[i];
       	rays.clear();
      }
      delete sep;
      num_lp_solved_++;
      OsiClpSolverInterface::resolve();
      // update soco_status_
      update_problem_status();
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
    num_lp_solved_++;
    OsiClpSolverInterface::resolve();
    // update problem status
    update_problem_status();
    // check problem status
    if (soco_status_!=OPTIMAL)
      break;
    delete sep;
    // chec if the new solution is feasible for cones
    sep = new Separate(cc_, getColSolution(), getNumCols(), options_);
    feasible = sep->is_feasible();
  }
  // update problem status
  update_problem_status();
  return soco_status_;
}

const ConicConstraints * ColaModel::get_conic_constraints() const{
  return cc_;
}

Options * ColaModel::options() const {
  return options_;
}


void ColaModel::setConicConstraints(ConicConstraints * cc) {
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
  const double * full_sol = getColSolution();
  for (int i=0; i<num_cones; ++i) {
    int cone_size = cc_->cone_size(i);
    const int * members = cc_->cone_members(i);
    par_sol = new double[cone_size];
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
  solve(false);
}

void ColaModel::resolve() {
  solve(true);
}

// returns problem status and updates status_
ProblemStatus ColaModel::update_problem_status() {
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

int ColaModel::getNumCones() const {
  int nc = cc_->num_cones();
  return nc;
}

void ColaModel::getConicConstraint(int index, OsiConeType & type,
				   int & numMembers,
				   int *& members) const {
  ConeType t = cc_->type(index);
  if (t==QUAD)
    type = OSI_QUAD;
  else
    type = OSI_RQUAD;
  numMembers = cc_->cone_size(index);
  //int nc = cc_->num_cones();
  members = new int[numMembers];
  const int * m = cc_->cone_members(index);
  std::copy(m, m+numMembers, members);
}

void ColaModel::removeConicConstraint(int index) {
  cc_->remove_cone(index);
  //num_cones_ = cc_->num_cones();
}

// first reduces all conic constraints to 3 dimentional second order conic
// constraints as described in ben-tal nemirovski, then solves the conic
// problem using liner approximations
ProblemStatus ColaModel::solve_reducing_cones(bool resolve) {
  int num_cones = cc_->num_cones();
  // if we have rotated cones bail out,
  for (int i=0; i<num_cones; ++i) {
    if (cc_->type(i)==OSI_RQUAD) {
      std::cerr << "Cola: This method is only for canonical cones."
		<< std::endl;
      std::cerr << "Cola: Terminating..." << std::endl;
      soco_status_ = ABANDONED;
      return ABANDONED;
    }
  }

  // create approximation
  // = first reduce all cones to 3 dimensional cones.
  // conic constraints obtained by reducing cone i
  ConicConstraints * reduced_cc_i = new ConicConstraints();
  // set of conic constraints we get by reducing all cones of original problem
  ConicConstraints * reduced_cc = new ConicConstraints();

  int num_var = getNumCols();
  for (int i=0; i<num_cones; ++i) {
    if (cc_->cone_size(i)<=3) {
      continue;
    }
    // reduce conic constraint i to smaller cones, save them in reduced_cc_i
    reduce_cone (cc_->cone_size(i), cc_->cone_members(i),
		 reduced_cc_i, num_var);
    // add new cones to problem
    int num_cones_i = reduced_cc_i->num_cones();
    for (int j=0; j<num_cones_i; ++j) {
      reduced_cc->add_cone(reduced_cc_i->cone_size(j),
			   reduced_cc_i->cone_members(j),
			   reduced_cc_i->type(j));
    }
    // reset reduced_cc_i
    reduced_cc_i->reset();
  }
  // print new cones of the problem
  //reduced_cc->dump_cones();
  // = add new variables to the model
  int diff = num_var - getNumCols();
  double infinity = getInfinity();
  for(int i=0; i<diff; ++i) {
    addCol(0, NULL, NULL, 0.0, infinity, 0.0);
  }
  // int num_rows = getNumRows();
  // int diff = num_var - getNumCols();
  // double * collb = new double[diff]();
  // double * colub = new double[diff]();
  // double * obj = new double[diff]();
  // double infinity = getInfinity();
  // std::fill(colub, colub+diff, infinity);
  // CoinPackedVectorBase ** cols = 0;
  // cols = new CoinPackedVectorBase*[diff];
  // int * inds = 0;
  // for (int i=0; i<diff; ++i) {
  //   cols[i] = new CoinPackedVector(num_rows, inds, 0.0);
  // }
  // addCols(diff, cols, collb, colub, obj);
  // std::cout << "Num Cols: " << getNumCols();
  // delete[] collb;
  // delete[] colub;
  // delete[] obj;
  //writeMps("reduced", "mps");
  // if reduced cone is not empty
  if (reduced_cc->num_cones()!=0) {
    setConicConstraints(reduced_cc);
  }
  int nOfCones = getNumCones();
  if (num_cuts_!=0) {
    delete[] num_cuts_;
  }
  if (num_supports_!=0) {
    delete num_supports_;
  }
  num_cuts_ = new int[nOfCones]();
  num_supports_ = new int[nOfCones]();
  solve(resolve);
  // set columns for new variables to 0
}

// solves problem using ben-tal nemirovski approximation
// v is approximation parameter
// resolve is false, then solve from scratch
// resolve is true, use Clp's resolve.
ProblemStatus ColaModel::solve_with_bn(bool resolve, int v) {
  // if we have rotated cones bail out,


  // create approximation
  // = first reduce all cones to 3 dimensional cones.
  // conic constraints obtained by reducing cone i
  ConicConstraints * reduced_cc_i = new ConicConstraints();
  // set of conic constraints we get by reducing all cones of original problem
  ConicConstraints * reduced_cc = new ConicConstraints();
  int num_cones = cc_->num_cones();
  int num_var = getNumCols();
  for (int i=0; i<num_cones; ++i) {
    if (cc_->cone_size(i)<=3) {
      continue;
    }
    // reduce conic constraint i to smaller cones, save them in reduced_cc_i
    reduce_cone (cc_->cone_size(i), cc_->cone_members(i),
		 reduced_cc_i, num_var);
    // add new cones to problem
    int num_cones_i = reduced_cc_i->num_cones();
    for (int j=0; j<num_cones_i; ++j) {
      reduced_cc->add_cone(reduced_cc_i->cone_size(j),
			   reduced_cc_i->cone_members(j),
			   reduced_cc_i->type(j));
    }
    // reset reduced_cc_i
    reduced_cc_i->reset();
  }
  // print new cones of the problem
  //reduced_cc->dump_cones();
  // = add new variables to the model
  int diff = num_var - getNumCols();
  double infinity = getInfinity();
  for(int i=0; i<diff; ++i) {
    addCol(0, NULL, NULL, 0.0, infinity, 0.0);
  }
  // int num_rows = getNumRows();
  // int diff = num_var - getNumCols();
  // double * collb = new double[diff]();
  // double * colub = new double[diff]();
  // double * obj = new double[diff]();
  // double infinity = getInfinity();
  // std::fill(colub, colub+diff, infinity);
  // CoinPackedVectorBase ** cols = 0;
  // cols = new CoinPackedVectorBase*[diff];
  // int * inds = 0;
  // for (int i=0; i<diff; ++i) {
  //   cols[i] = new CoinPackedVector(num_rows, inds, 0.0);
  // }
  // addCols(diff, cols, collb, colub, obj);
  // std::cout << "Num Cols: " << getNumCols();
  // delete[] collb;
  // delete[] colub;
  // delete[] obj;
  //writeMps("reduced", "mps");
  // if reduced cone is not empty
  if (reduced_cc->num_cones()!=0) {
    setConicConstraints(reduced_cc);
  }
  int nOfCones = getNumCones();
  if (num_cuts_!=0) {
    delete[] num_cuts_;
  }
  if (num_supports_!=0) {
    delete num_supports_;
  }
  num_cuts_ = new int[nOfCones]();
  num_supports_ = new int[nOfCones]();
  solve(resolve);
  // set columns for new variables to 0
}

// this function assumes cone is in canonical form.
void ColaModel::reduce_cone(int size, const int * members,
			    ConicConstraints * reduced_cc_i, int & num_var) {
  // k is
  if (size<=3) {
    // add current cone and return
    reduced_cc_i->add_cone(size, members, QUAD);
    return;
  }
  int k = size/2;
  int * nm = new int[3];
  for (int i=0; i<k-1; ++i) {
    nm[0] = num_var+i;
    nm[1] = members[2*i+1];
    nm[2] = members[2*i+2];
    reduced_cc_i->add_cone(3, nm, QUAD);
  }
  // check if the last cone has 3 or 2 members
  if (2*k==size) {
    nm[0] = num_var+k-1;
    nm[1] = members[size-1];
    reduced_cc_i->add_cone(2, nm, QUAD);
  }
  else {
    nm[0] = num_var+k-1;
    nm[1] = members[size-2];
    nm[2] = members[size-1];
    reduced_cc_i->add_cone(3, nm, QUAD);
  }
  // reduce cone of new variables
  int new_cone[k+1];
  new_cone[0] = members[0];
  for (int i=1; i<k+1; ++i) {
    new_cone[i] = num_var+i-1;
  }
  num_var = num_var+k;
  // recursive call for reducing the new cone.
  reduce_cone(k+1, new_cone, reduced_cc_i, num_var);
}

int ColaModel::num_lp_solved() const {
  return num_lp_solved_;
}

ProblemStatus ColaModel::problem_status() const {
  return soco_status_;
}
