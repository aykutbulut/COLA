#include "ColaConfig.h"
#include "ColaModel.hpp"
#include "Cone.hpp"
#include "LorentzCone.hpp"
#include "ScaledCone.hpp"

#include "Separate.hpp"

#include "CoinMpsIO.hpp"

#include <OsiClpSolverInterface.hpp>
#include <CoinPackedMatrix.hpp>
#include <CoinPackedVector.hpp>

#include <cstring>
#include <iomanip>
#include <numeric>
#include <cmath>
#include <string>
#include <iostream>

#ifndef ipfint
//typedef ipfint FORTRAN_INTEGER_TYPE ;
typedef int ipfint;
typedef const int cipfint;
#endif

// using simple lapack interface
extern "C"
{
  /** LAPACK Fortran subroutine DPOSV. */
  void F77_FUNC(dposv, DPOSV)(char * uplo, ipfint * n, ipfint * nrhs,
                              double * A, ipfint * lda, double * B,
                              ipfint * ldb, ipfint * info);
}

ColaModel::ColaModel() : OsiClpSolverInterface() {
  options_ = new Options();
  total_num_supports_ = 0;
  total_num_cuts_ = 0;
  num_supports_ = 0;
  num_cuts_ = 0;
  num_lp_solved_ = 0;
  imp_solution_ = 0;
  cones_.clear();
  // this is the default beavior, user can change this using options
  setHintParam(OsiDoReducePrint,true,OsiHintTry);
  // for unboundedness directions set option
  OsiClpSolverInterface::getModelPtr()->setMoreSpecialOptions(0);
}

ColaModel::ColaModel(char * data_file) : OsiClpSolverInterface() {
  options_ = new Options();
  total_num_supports_ = 0;
  total_num_cuts_ = 0;
  num_supports_ = 0;
  num_cuts_ = 0;
  num_lp_solved_ = 0;
  imp_solution_ = 0;
  // this is the default beavior, user can change this using options
  setHintParam(OsiDoReducePrint,true,OsiHintTry);
  cones_.clear();
  OsiConicSolverInterface::readMps(data_file);
  OsiClpSolverInterface::getModelPtr()->setMoreSpecialOptions(0);
}

// copy constructor
ColaModel::ColaModel(ColaModel const & other): OsiClpSolverInterface(other) {
  // copy conic constraints
  std::vector<Cone*> other_cones = other.get_conic_constraints();
  std::vector<Cone*>::const_iterator it;
  for (it=other_cones.begin(); it!=other_cones.end(); ++it) {
    cones_.push_back((*it)->clone());
  }
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
  imp_solution_ = 0;
}

// copy assignment operator
ColaModel & ColaModel::operator=(ColaModel const & rhs) {
  // copy conic constraints
  std::vector<Cone*> rhs_cones = rhs.get_conic_constraints();
  std::vector<Cone*>::const_iterator it;
  for (it=rhs_cones.begin(); it!=rhs_cones.end(); ++it) {
    cones_.push_back((*it)->clone());
  }
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
  imp_solution_ = 0;
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
  // free conic constraints
  // std::vector<Cone*>::iterator it;
  // for (it=cones_.begin(); it!=cones_.end(); it++) {
  //   delete *it;
  // }
  cones_.clear();
  // free rest
  if(num_cuts_) {
    delete[] num_cuts_;
  }
  if (num_supports_) {
    delete[] num_supports_;
  }
  if (imp_solution_) {
    delete[] imp_solution_;
  }
}


void ColaModel::addConicConstraint(OsiLorentzConeType type,
				  int numMembers,
				  int const * members) {
  ConeType t;
  if (type==OSI_QUAD)
    t = LORENTZ;
  else
    t = RLORENTZ;
  Cone * c = new LorentzCone(t, numMembers, members);
  cones_.push_back(c);
}


void ColaModel::addConicConstraint(CoinPackedMatrix const * A,
				   CoinPackedVector const * b,
				   CoinPackedVector const * d, double h) {
  Cone * c = new ScaledCone(A, b, d, h);
  cones_.push_back(c);
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
  if (num_cuts_)
    delete[] num_cuts_;
  if (num_supports_)
    delete[] num_supports_;
  num_cuts_ = new int[nOfCones]();
  num_supports_ = new int[nOfCones]();
}

void ColaModel::print_stats() const {
  // how many separ problem solved
  // separation for each cone
  int num_cones = cones_.size();
  std::cout << "Total Unboundedness supports: " << total_num_supports_ << std::endl;
  if (total_num_supports_!=0) {
    for (int i=0; i<num_cones; ++i) {
      std::cout << "     Supports for cone " << std::setw(5) << i
		<< std::setw(6) << num_supports_[i] << std::endl;
    }
  }
  std::cout << "Total Seperating hyperplanes: " << total_num_cuts_ << std::endl;
  if (total_num_cuts_) {
    for (int i=0; i<num_cones; ++i) {
      std::cout << "     Seperating hyperplane for cone " << std::setw(5) << i
		<< std::setw(6) << num_cuts_[i] << std::endl;
    }
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
  // todo(aykut) now we should check if the cone related data is initialized
  // properly. This is a problem when the user does not read problem from mps
  // file but builds the model herself using cola as a library.
  if (num_cuts_==0)
    initialSolve();

  // solve linear problem
  bool feasible = false;
  Separate * sep;
  // add nonnegativity of leading variables
  std::vector<Cone*>::const_iterator it;
  for (it=cones_.begin(); it!=cones_.end(); ++it) {
    (*it)->relax(*this);
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
    if (options_->get_int_option(LOG_LEVEL)>0) {
      std::cout << "Cola: Problem without conic constraints is unbounded. Adding"
        " supporting hyperplanes to resolve..."
                << std::endl;
    }
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
      // std::cout << "Unboundedness direction is " << std::endl  << "[";
      // for (int i=0; i<getNumCols(); ++i) {
      // 	std::cout << std::setw(10) << vec[i] << "; ";
      // }
      // std::cout << "]" << std::endl;
      // END OF DIRECTION PRINT
      sep = new Separate(cones_, vec, getNumCols(), options_);
      // returns true if direction is feasible for all cones.
      feasible = sep->is_feasible();
      if (feasible) {
	// primal ray is feasible for all cone constraints,
	// problem is unbounded
	delete sep;
	soco_status_ = DUAL_INFEASIBLE;
	return DUAL_INFEASIBLE;
      }
      else {
	// Add all the cuts generated
	std::vector<CoinPackedVector*> cut = sep->cuts();
	std::vector<double> rhs = sep->rhs();
	std::vector<int> gen_cone = sep->generating_cone();
	int num_cuts = cut.size();
	for (int i=0; i<num_cuts; ++i) {
	  addRow(*cut[i], -getInfinity(), rhs[i]);
	  // print cut
	  // end of cut print
	  total_num_supports_++;
	  num_supports_[gen_cone[i]]++;
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
  if (soco_status_!=OPTIMAL) {
    // it means it is not optimal after resolving unbounded directions
    std::cout << "Cola: Problem status is not optimal after adding "
      "supporting hyperplanes."
	      << std::endl;
    std::cout << "Cola: Terminating..." << std::endl;
    return soco_status_;
  }
  // when we reached here, it means the problem is not unbounded anymore.
  // todo(aykut) pick the cone and all those logs and stuff. it can all go to separate.
  sep = new Separate(cones_, getColSolution(), getNumCols(), options_);
  // returns true if given point is feasible for all cones.
  feasible = sep->is_feasible();
  while(!feasible) {
    // number of cuts generated
    // std::cout << "ColaModel: " << sep->cuts()->size() << " cuts generated." << std::endl;
    // add all the cuts generated
    // Add all the cuts generated
    std::vector<CoinPackedVector*> cut = sep->cuts();
    std::vector<double> rhs = sep->rhs();
    std::vector<int> gen_cone = sep->generating_cone();
    int num_cuts = cut.size();
    for (int i=0; i<num_cuts; ++i) {
      addRow(*cut[i], -getInfinity(), rhs[i]);
      // print cut
      // end of cut print
      total_num_cuts_++;
      num_cuts_[gen_cone[i]]++;
    }
    // if (num_lp_solved_%2==0) {
    //   clean_redundant_constraints();
    // }
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
    sep = new Separate(cones_, getColSolution(), getNumCols(), options_);
    feasible = sep->is_feasible();
  }
  delete sep;
  // update problem status
  update_problem_status();
  // clean redundant constraints
  // clean_redundant_constraints();
  return soco_status_;
}

std::vector<Cone*> ColaModel::get_conic_constraints() const{
  return cones_;
}

Options * ColaModel::options() const {
  return options_;
}


void ColaModel::setConicConstraints(std::vector<Cone*> cones) {
  cones_ = cones;
}

void ColaModel::setOptions(Options * opt) {
  options_ = opt;
}

void ColaModel::report_feasibility() const {
  std::cout << "Conic Constraints feasibility report" << std::endl;
  std::cout << std::setw(5) << std::left << "Cone";
  // todo(aykut) this is not true all the time, what if cone is rotated.
  std::cout << std::setw(20) << std::left << "x1^2-sum x_i^2"
	    << std::setw(20) << std::left << "x1-||x_{2:n}||"
	    << std::endl;
  for (int i=0; i<cones_.size(); ++i) {
    if (cones_[i]->type()==LORENTZ) {
      std::cout << std::setw(5) << std::left << i
		<< std::setw(20) << std::left << "-"
		<< std::setw(20) << std::left << cones_[i]->feasibility(getColSolution())
		<< std::endl;
    }
    else if (cones_[i]->type()==RLORENTZ) {
      std::cout << std::setw(5) << std::left << i
		<< std::setw(20) << std::left << cones_[i]->feasibility(getColSolution())
		<< std::setw(20) << std::left << "-"
		<< std::endl;
    }
  }
}

void ColaModel::initialSolve() {
  int nOfCones = getNumCones();
  if (num_cuts_==0) {
    num_cuts_ = new int[nOfCones]();
  }
  if (num_supports_==0) {
    num_supports_ = new int[nOfCones]();
  }
  solve(false);
  update_problem_status();
}

void ColaModel::resolve() {
  solve(true);
  update_problem_status();
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
  int nc = cones_.size();
  std::copy(num_cuts, num_cuts+nc, num_cuts_);
}

void ColaModel::set_num_supports(const int * num_supports) {
  int nc = cones_.size();
  std::copy(num_supports, num_supports+nc, num_supports_);
}

void ColaModel::set_total_num_cuts(int tnc) {
  total_num_cuts_ = tnc;
}

void ColaModel::set_total_num_supports(int tns) {
  total_num_supports_ = tns;
}

int ColaModel::getNumCones() const {
  int nc = cones_.size();
  return nc;
}

int ColaModel::getConeSize(int i) const {
  return cones_.size();
}

OsiConeType ColaModel::getConeType(int i) const {
  ConeType type;
  type = cones_[i]->type();
  if (type==SCALED)
    return OSI_SCALED;
  else
    return OSI_LORENTZ;
}

void ColaModel::getConeSize(int * size) const {
  int num_cones = cones_.size();
  for (int i=0; i<num_cones; ++i) {
    size[i] = cones_[i]->size();
  }
}

void ColaModel::getConeType(OsiConeType * type) const {
  int num_cones = cones_.size();
  ConeType t;
  for (int i=0; i<num_cones; ++i) {
    t = cones_[i]->type();
    if (t==SCALED)
      type[i] = OSI_SCALED;
    else
      type[i] = OSI_LORENTZ;
  }
}

void ColaModel::getConeType(OsiLorentzConeType * type) const {
  int num_cones = cones_.size();
  ConeType t;
  for (int i=0; i<num_cones; ++i) {
    t = cones_[i]->type();
    if (t==LORENTZ)
      type[i] = OSI_QUAD;
    else if (t==RLORENTZ)
      type[i] = OSI_RQUAD;
    else {
      std::cerr << "Cone is not a lorentz cone!" << std::endl;
      throw std::exception();
    }
  }
}

void ColaModel::getConicConstraint(int index, OsiLorentzConeType & type,
				   int & numMembers,
				   int *& members) const {
  ConeType t = cones_[index]->type();
  if (t==SCALED) {
    std::cerr << "this function is for Lorentz cones!" << std::endl;
    throw std::exception();
  }
  if (t==LORENTZ)
    type = OSI_QUAD;
  else if (t==RLORENTZ)
    type = OSI_RQUAD;
  numMembers = cones_[index]->size();
  members = new int[numMembers];
  LorentzCone * c = dynamic_cast<LorentzCone*>(cones_[index]);
  const int * m = c->members();
  std::copy(m, m+numMembers, members);
}

void ColaModel::removeConicConstraint(int index) {
  cones_.erase(cones_.begin()+index);
}

void ColaModel::modifyConicConstraint(int index, OsiLorentzConeType type,
				      int numMembers,
				      int const * members) {
  // free cone data first
  delete cones_[index];
  // insert new cone
  ConeType t;
  if (type==OSI_QUAD)
    t = LORENTZ;
  else
    t = RLORENTZ;
  cones_[index] = new LorentzCone(t, numMembers, members);
}

// first reduces all conic constraints to 3 dimentional second order conic
// constraints as described in ben-tal nemirovski, then solves the conic
// problem using liner approximations
ProblemStatus ColaModel::solve_reducing_cones(bool resolve) {
  int num_cones = cones_.size();
  // if we have Scaled cones bail out,
  std::vector<Cone*>::const_iterator it;
  for (it=cones_.begin(); it!=cones_.end(); it++) {
    if ((*it)->type()==SCALED) {
      std::cerr << "Cola: This method is only for Lorentz cones."
		<< std::endl;
      std::cerr << "Cola: Terminating..." << std::endl;
      soco_status_ = ABANDONED;
      return ABANDONED;
    }
  }
  // if we have rotated cones bail out,
  for (int i=0; i<num_cones; ++i) {
    if (cones_[i]->type()==SCALED or cones_[i]->type()==RLORENTZ) {
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
  std::vector<Cone*> reduced_cones_i;
  // set of conic constraints we get by reducing all cones of original problem
  std::vector<Cone*>  reduced_cones;
  int num_var = getNumCols();
  for (int i=0; i<num_cones; ++i) {
    LorentzCone * c = dynamic_cast<LorentzCone*>(cones_[i]);
    if (c->size()<=3) {
      continue;
    }
    // reduce conic constraint i to smaller cones, save them in reduced_cc_i
    reduce_cone (c->size(), c->members(), reduced_cones_i, num_var);
    // add new cones to problem
    int num_cones_i = reduced_cones_i.size();
    for (int i=0; i<num_cones_i; ++i) {
      reduced_cones.push_back(reduced_cones_i[i]);
    }
    // copy all from reduced_cones_i vector to reduced_cones
    //    std::copy(reduced_cones_i, reduced_cones_i+reduced_cones_i.size(),
    //	      reduced_cones+reduced_cones.size());
    // reset reduced_cc_i
    reduced_cones_i.clear();;
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
  if (!reduced_cones.empty()) {
    setConicConstraints(reduced_cones);
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
  std::vector<Cone*> reduced_cones_i;
  // set of conic constraints we get by reducing all cones of original problem
  std::vector<Cone*> reduced_cones;
  int num_cones = cones_.size();
  int num_var = getNumCols();
  for (int i=0; i<num_cones; ++i) {
    if (cones_[i]->size()<=3) {
      continue;
    }
    // reduce conic constraint i to smaller cones, save them in reduced_cc_i
    LorentzCone * c = dynamic_cast<LorentzCone*>(cones_[i]);
    reduce_cone (c->size(), c->members(), reduced_cones_i, num_var);
    // add new cones to problem
    int num_cones_i = reduced_cones_i.size();
    for (int i=0; i<num_cones_i; ++i) {
      reduced_cones.push_back(reduced_cones_i[i]);
    }
    //    std::copy(reduced_cones_i, reduced_cones_i+reduced_cones_i.size(),
    //	      reduced_cones+reduced_cones.size());
    // reset reduced_cc_i
    reduced_cones_i.clear();
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
  if (!reduced_cones.empty()) {
    setConicConstraints(reduced_cones);
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
			    std::vector<Cone*> & reduced_cones_i, int & num_var) {
  // k is
  if (size<=3) {
    // add current cone and return
    Cone * c = new LorentzCone(LORENTZ, size, members);
    reduced_cones_i.push_back(c);
    return;
  }
  int k = size/2;
  int * nm = new int[3];
  for (int i=0; i<k-1; ++i) {
    nm[0] = num_var+i;
    nm[1] = members[2*i+1];
    nm[2] = members[2*i+2];
    Cone * c = new LorentzCone(LORENTZ, 3, nm);
    reduced_cones_i.push_back(c);
  }
  // check if the last cone has 3 or 2 members
  if (2*k==size) {
    nm[0] = num_var+k-1;
    nm[1] = members[size-1];
    Cone * c = new LorentzCone(LORENTZ, 2, nm);
    reduced_cones_i.push_back(c);
  }
  else {
    nm[0] = num_var+k-1;
    nm[1] = members[size-2];
    nm[2] = members[size-1];
    Cone * c = new LorentzCone(LORENTZ, 3, nm);
    reduced_cones_i.push_back(c);
  }
  // reduce cone of new variables
  int new_cone[k+1];
  new_cone[0] = members[0];
  for (int i=1; i<k+1; ++i) {
    new_cone[i] = num_var+i-1;
  }
  num_var = num_var+k;
  // recursive call for reducing the new cone.
  reduce_cone(k+1, new_cone, reduced_cones_i, num_var);
}

int ColaModel::num_lp_solved() const {
  return num_lp_solved_;
}

ProblemStatus ColaModel::problem_status() const {
  return soco_status_;
}

void ColaModel::clean_redundant_constraints() {
  int num_rows = getNumRows();
  int num_cols = getNumCols();
  int cuts_added = total_num_cuts_ + total_num_supports_;
  int orig_num_rows = num_rows - total_num_cuts_;
  std::vector<int> rows_to_delete;
  int * cstat = new int[num_cols];
  int * rstat = new int[num_rows];
  // getBasisStatus(cstat, rstat);
  double const * row_act = getRowActivity();
  for (int i=orig_num_rows; i<num_rows-20; ++i) {
    if (rstat[i]==1)
      rows_to_delete.push_back(i);
    // delete rows that are irrelevant
    // check the activity, if it is too redundant less than 1e-5,
    // then remove it.
    // if (abs(row_act[i]) > 0.001*options_->get_dbl_option(TOL)) {
    //   rows_to_delete.push_back(i);
    // }
  }
  int num_to_delete = rows_to_delete.size();
  int * row_indices = new int[num_to_delete];
  std::copy(rows_to_delete.begin(), rows_to_delete.end(), row_indices);
  if (num_to_delete) {
    OsiClpSolverInterface::deleteRows(num_to_delete, row_indices);
    // update total number of cuts
    total_num_cuts_ -=  num_to_delete;
  }
  delete[] cstat;
  delete[] rstat;
  delete[] row_indices;
  if (num_to_delete)
    std::cout << num_to_delete << " redundant constraints removed."
	      << std::endl;
}

void ColaModel::writeMps(const char * filename,
			 const char * extension,
			 double objSense) const {
  std::cerr << "Writing mps files are not implemented yet." << std::endl;
  OsiClpSolverInterface::writeMps(filename, extension, objSense);
  //throw std::exception();
}

ProblemStatus ColaModel::solve_numeric() {
  solve(false);
  // 1. determine original binding constraints
  std::vector<int> binding;
  int num_rows = getNumRows();
  int orig_num_rows = num_rows - total_num_cuts_;
  int num_cols = getNumCols();
  int * cstat = new int[num_cols];
  int * rstat = new int[num_rows];
  getBasisStatus(cstat, rstat);
  for (int i=0; i<orig_num_rows; ++i) {
    if (rstat[i]==2 or rstat[i]==3) {
      binding.push_back(i);
    }
  }
  int num_binding = binding.size();
  delete[] cstat;
  delete[] rstat;
  // 2. compute AA^T
  // A is sparse and instance of CoinPackedMatrix
  // 2.1 get A
  CoinPackedMatrix const * mat = getMatrixByRow();
  CoinPackedMatrix * A = new CoinPackedMatrix(false, 500, 500);
  std::vector<int>::const_iterator it;
  for (it=binding.begin(); it!=binding.end(); ++it) {
    int first = mat->getVectorFirst(*it);
    int last = mat->getVectorLast(*it);
    int num_elem  = last-first;
    int const * mat_cols = mat->getIndices() + first;
    double const * mat_value = mat->getElements() + first;
    A->appendRow(CoinPackedVector(num_elem, mat_cols, mat_value));
  }
  // 2.2 Compute AAt
  CoinPackedMatrix * AAt = new CoinPackedMatrix();
  double * Aa_i = new double[num_binding];
  int * Aa_i_cols = new int[num_binding];
  double * Aa_i_vals = new double[num_binding];
  for (it=binding.begin(); it!=binding.end(); ++it) {
    // A times row i of A
    int first = mat->getVectorFirst(*it);
    int last = mat->getVectorLast(*it);
    int num_elem  = last-first;
    int const * mat_cols = mat->getIndices() + first;
    double const * mat_value = mat->getElements() + first;
    A->times(CoinPackedVector(num_elem, mat_cols, mat_value), Aa_i);
    // sparsify and insert Aa_i
    int Aa_i_size = 0;
    for (int i=0; i<num_binding; ++i) {
      if (Aa_i[i]!=0.0) {
	Aa_i_cols[Aa_i_size] = i;
	Aa_i_vals[Aa_i_size] = Aa_i[i];
	Aa_i_size++;
      }
    }
    AAt->appendCol(CoinPackedVector(Aa_i_size, Aa_i_cols, Aa_i_vals));
  }
  delete[] Aa_i;
  delete[] Aa_i_cols;
  delete[] Aa_i_vals;
  //AAt->dumpMatrix();
  // 3. compute Ac
  double * Ac = new double[num_binding];
  double const * obj = getObjCoefficients();
  int obj_size;
  int * obj_cols = new int[num_cols];
  double * obj_vals = new double[num_cols];
  for (int i=0; i<num_cols; ++i) {
    if (obj[i]!=0) {
      obj_cols[obj_size] = i;
      obj_vals[obj_size] = obj[i];
      obj_size++;
    }
  }
  A->times(CoinPackedVector(obj_size, obj_cols, obj_vals), Ac);
  delete[] obj_cols;
  delete[] obj_vals;
  // 4. compute (b-Ac)
  double * b = new double[num_binding];
  int k=0;
  for (it=binding.begin(); it!=binding.end(); ++it) {
    b[k] = getRightHandSide()[*it];
    k++;
  }
  double * b_Ac = new double[num_binding];
  for (int i=0; i<num_binding; ++i) {
    b_Ac[i] = b[i] - Ac[i];
  }
  // 5. solve AA^Ty=(b-Ac)
  // 5.1 get AAt in lower triangular format
  double ** AAt_dense = new double*[num_binding];
  for (int i=0; i<num_binding; ++i) {
    AAt_dense[i] = new double[num_binding]();
  }
  int const * AAt_cols = AAt->getIndices();
  double const * AAt_value = AAt->getElements();
  for (int i=0; i<num_binding; ++i) {
    // get row i
    int first = AAt->getVectorFirst(i);
    int last = AAt->getVectorLast(i);
    //int num_elem  = last-first;
    for (int j=first; j<last; ++j) {
      AAt_dense[i][AAt_cols[j]] = AAt_value[j];
    }
  }
  // 5.2 call lapack routine to solve the system
  double * y = new double[num_binding];
  lapack_solve(AAt_dense, b_Ac, y, num_binding);
  // 6. compute d <- c+A^Ty
  // in other words x <- c + A'(AA')^{-1}b - A'(AA')^{-1}Ac when
  // we insert for y.
  // 6.1 compute Aty
  double * Aty = new double[num_cols];
  A->transposeTimes(y, Aty);
  // 6.2 compute d <- c+A^Ty
  double * dir = new double[num_cols];
  double const * cur_sol = getColSolution();
  for (int i=0; i<num_cols; ++i) {
    dir[i] = -obj[i] - Aty[i];
  }
  // 7. go along direction until we hit boundry, ie. compute step_size
  double step_size = 0.0;
  int num_cones = getNumCones();
  imp_solution_ = new double[num_cols];
  std::copy(cur_sol, cur_sol+num_cols, imp_solution_);
  double * par_sol;
  double * par_dir;
  for (int i=0; i<num_cones; ++i) {
    int cone_size = cones_[i]->size();
    LorentzCone * con = dynamic_cast<LorentzCone*>(cones_[i]);
    int const * members = con->members();
    par_sol = new double[cone_size];
    par_dir = new double[cone_size];
    // 7.1 compute step size
    // 7.2 compute a in ax^2+bx+c=0
    // 7.3 compute b in ax^2+bx+c=0
    // 7.4 compute c in ax^2+bx+c=0
    for (int j=0; j<con->size(); j++) {
      par_sol[j] = cur_sol[members[j]];
      par_dir[j] = dir[members[j]];
    }
    double feas_i = par_sol[0]*par_sol[0]-std::inner_product(par_sol+1, par_sol+cone_size, par_sol+1, 0.0);
    if (feas_i > options_->get_dbl_option(TOL)) {
      delete[] par_sol;
      delete[] par_dir;
      continue;
    }
    double a = par_dir[0]*par_dir[0];
    double b = par_sol[0]*par_dir[0];
    double c = par_sol[0]*par_sol[0];
    a = a - std::inner_product(par_dir+1, par_dir+con->size(), par_dir+1, 0.0);
    b = b - std::inner_product(par_sol+1, par_sol+con->size(), par_dir+1, 0.0);
    b = 2.0*b;
    c = c - std::inner_product(par_sol+1, par_sol+con->size(), par_sol+1, 0.0);
    // 7.5 compute step size
    double alpha1 = (-b+sqrt(b*b-4.0*a*c))/(2.0*a);
    double alpha2 = (-b-sqrt(b*b-4.0*a*c))/(2.0*a);
    double alpha=alpha1;
    std::cout << "alpha1 " << alpha1 << std::endl;
    std::cout << "alpha2 " << alpha2 << std::endl;
    // if (alpha2<alpha1)
    //   alpha=alpha2;
    for (int j=0; j<con->size(); j++) {
      imp_solution_[members[j]] = cur_sol[members[j]] + alpha*par_dir[j];
    }
    delete[] par_sol;
    delete[] par_dir;
    // get related portion of direction
    // get related portion of solution
  }
  delete A;
  delete AAt;
  delete[] Ac;
  delete[] b;
  delete[] b_Ac;
  for (int i=0; i<num_binding; ++i) {
    delete[] AAt_dense[i];
  }
  delete[] AAt_dense;
  delete[] y;
  delete[] Aty;
  delete[] dir;
  // remove all the cuts and add a cut using the improved solution
  int * indices = new int[total_num_cuts_];
  for (int i=0; i<total_num_cuts_; ++i) {
    indices[i] = orig_num_rows+i;
  }
  deleteRows(total_num_cuts_, indices);
  delete[] indices;
  total_num_cuts_ = 0;
  num_lp_solved_ = 0;
  std::fill(num_cuts_, num_cuts_+num_cones, 0);
  // separate from imp_solution
  Separate * sep = new Separate(cones_, imp_solution_, getNumCols(), options_);
  // returns true if given point is feasible for all cones.
  int feasible = sep->is_feasible();
  while(!feasible) {
    // number of cuts generated
    // std::cout << "ColaModel: " << sep->cuts()->size() << " cuts generated." << std::endl;
    // add all the cuts generated
	// Add all the cuts generated
    std::vector<CoinPackedVector*> cut = sep->cuts();
    std::vector<double> rhs = sep->rhs();
    std::vector<int> gen_cone = sep->generating_cone();
    int num_cuts = cut.size();
    for (int i=0; i<num_cuts; ++i) {
      addRow(*cut[i], -getInfinity(), rhs[i]);
      // print cut
      // end of cut print
      total_num_cuts_++;
      num_cuts_[gen_cone[i]]++;
    }
    // if (num_lp_solved_%2==0) {
    //   clean_redundant_constraints();
    // }
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
    sep = new Separate(cones_, getColSolution(), getNumCols(), options_);
    feasible = sep->is_feasible();
  }
  delete sep;
  // update problem status
  update_problem_status();
  // report feasibility of imp_solution_
  double lhs = 0.0;
  double lhs_real = 0.0;
  std::cout << "Conic Constraints feasibility report" << std::endl;
  std::cout << std::setw(5) << std::left << "Cone";
  // todo(aykut) this is not true all the time, what if cone is rotated.
  std::cout << std::setw(20) << std::left << "x1^2-sum x_i^2"
	    << std::setw(20) << std::left << "x1-||x_{2:n}||"
	    << std::endl;
  const double * full_sol = imp_solution_;
  //double * par_sol;
  for (int i=0; i<num_cones; ++i) {
    int cone_size = cones_[i]->size();
    LorentzCone * con = dynamic_cast<LorentzCone*>(cones_[i]);
    const int * members = con->members();
    par_sol = new double[cone_size];
    for (int j=0; j<cone_size; ++j) {
      par_sol[j] = full_sol[members[j]];
    }
    if (con->type()==LORENTZ) {
      lhs = par_sol[0]*par_sol[0]
	- std::inner_product(par_sol+1, par_sol+cone_size, par_sol+1, 0.0);
      lhs_real = par_sol[0]
	-sqrt(std::inner_product(par_sol+1, par_sol+cone_size, par_sol+1, 0.0));
    }
    else if (con->type()==RLORENTZ) {
      lhs = 2.0*par_sol[0]*par_sol[1]
	- std::inner_product(par_sol+2, par_sol+cone_size, par_sol+2, 0.0);
      lhs_real = lhs;
    }
    std::cout << std::setw(5) << std::left << i
              << std::setw(20) << std::left << lhs
              << std::setw(20) << std::left << lhs_real
              << std::endl;
    delete[] par_sol;
  }
  return soco_status_;
}

// solves system of symmertic positive definite Ax=b and store the solution
// in x.
void ColaModel::lapack_solve(double ** A, double * b, double * x, int dim) {
  char uplo = 'L';
  int num_rhs = 1;
  int info;
  // copy b to x
  std::copy(b, b+dim, x);
  int size = (dim*dim+dim)/2;
  double * A_lower = new double[size];
  for (int i=0; i<dim; ++i) {
    std::copy(A[i], A[i]+i+1, A_lower+(i*i+i)/2);
  }
  F77_FUNC (dposv, DPOSV) (&uplo, &dim, &num_rhs, A_lower, &dim, x, &dim, &info);
  if (info!=0) {
    std::cerr << "Lapack dposv function failed." << std::endl;
    throw std::exception();
  }
  delete[] A_lower;
}
