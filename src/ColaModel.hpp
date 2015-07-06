#ifndef COLA_MODEL_H
#define COLA_MODEL_H

//#include <OsiSolverInterface.hpp>
#include <OsiClpSolverInterface.hpp>
#include <OsiConicSolverInterface.hpp>
#include "Cone.hpp"
#include "Options.hpp"

typedef enum {
  ABANDONED=0,
  OPTIMAL,
  PRIMAL_INFEASIBLE,
  DUAL_INFEASIBLE,
  PRIMAL_OBJECTIVE_LIMIT_REACHED,
  DUAL_OBJECTIVE_LIMIT_REACHED,
  ITERATION_LIMIT_REACHED,
} ProblemStatus;

// TODO(aykut) implement getIterationCount(). It should give number of simplex
// iteration count. Right now it gives the number of iteration count of
// the last OsiClpSolverInterface::resolve call.

// TODO(aykut) how relevant the dual solution is? That is used by symphony.
//


// does the following make sense?
// todo(aykut) when cola solve method is invoked it should create its own copy
// of problem data and work on the copy. This is so since Cola modifies the
// problem by adding approximating hyperplanes.

class ColaModel: virtual public OsiConicSolverInterface,
		 public OsiClpSolverInterface {
  // data members
  std::vector<Cone*> cones_;
  Options * options_;
  ProblemStatus soco_status_;
  // number of times we solve lp relaxation problem
  int num_lp_solved_;
  // number of cuts generated for each cone
  int * num_cuts_;
  // total number of cuts
  int total_num_cuts_;
  // number of supports generated for each cone
  int * num_supports_;
  // total number of supports
  int total_num_supports_;
  double * imp_solution_;
  // PRIVATE FUNCTIONS
  // reduce conic constraint given by size and members to smaller cones, save
  // them in reduced_cone
  // num_var will store the new number of variables after reduction
  void reduce_cone(int size, const int * members,
		   std::vector<Cone*> & reduced_cone, int & num_var);
protected:
  // get cut and support statistics
  const int * num_cuts() const;
  int num_cuts(int i) const;
  const int * num_supports() const;
  int num_supports(int i) const;
  int total_num_cuts() const;
  int total_num_supports() const;
  // set cut and support statistics
  void set_num_cuts(const int * num_cuts);
  void set_num_supports(const int * num_supports);
  void set_total_num_cuts(int tnc);
  void set_total_num_supports(int tns);
  ProblemStatus solve(bool resolve);
  // solve symmetric positive definite linear system Ax=b
  void lapack_solve(double ** A, double * b, double * x, int dim);
public:
  // default constructor
  ColaModel();
  // constructor that reads mps file
  ColaModel(char * data_file);
  // copy constructor
  ColaModel(const ColaModel & other);
  // copy assignment operator
  ColaModel & operator=(const ColaModel & rhs);
  // destructor
  virtual ~ColaModel();
  // update and return problem status
  ProblemStatus update_problem_status();
  // return problem status;
  ProblemStatus problem_status() const;
  // functions
  //void read(const char * data_file);
  // following three is useful when we clone
  void setConicConstraints(std::vector<Cone*> cones);
  std::vector<Cone*> get_conic_constraints() const;
  void setOptions(Options * opt);
  Options * options() const;
  void print_stats() const;
  void print_solution() const;
  void report_feasibility() const;
  int num_lp_solved() const;
  // VIRTUAL FUNCTIONS
  // get conic constraints
  virtual void getConicConstraint(int index, OsiLorentzConeType & type,
				  int & numMembers,
				  int *& members) const;
  // add conic constraints
  // add conic constraint in lorentz cone form
  virtual void addConicConstraint(OsiLorentzConeType type,
				  int numMembers,
				  int const * members);
  // add conic constraint in |Ax-b| <= dx-h form
  virtual void addConicConstraint(CoinPackedMatrix const * A,
				  CoinPackedVector const * b,
				  CoinPackedVector const * d, double h);
  virtual void removeConicConstraint(int index);
  virtual void modifyConicConstraint(int index, OsiLorentzConeType type,
				     int numMembers,
				     int const * members);
  //ProblemStatus solve();
  virtual int getNumCones() const;
  virtual int getConeSize(int i) const;
  virtual void getConeSize(int * size) const;
  virtual OsiConeType getConeType(int i) const;
  virtual void getConeType(OsiConeType * type) const;
  virtual void getConeType(OsiLorentzConeType * type) const;
  virtual int readMps(const char * filename, const char * extension="mps");
  virtual void initialSolve();
  virtual void resolve();
  virtual void solveFromHotStart();
  virtual OsiConicSolverInterface * clone (bool copyData=true) const;
  // END OF VIRTUAL FUNCTIONS
  // solves problem using ben-tal nemirovski approximation
  // v is approximation paramete
  ProblemStatus solve_reducing_cones(bool resolve);
  ProblemStatus solve_with_bn(bool resolve, int v);
  ProblemStatus solve_numeric();
  void clean_redundant_constraints();
  // use conic solver interface's writeMps method
  virtual void writeMps(const char * filename,
			const char * extension = "mps",
			double objSense=0.0) const;
};

#endif
