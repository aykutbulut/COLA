#ifndef COLA_MODEL_H
#define COLA_MODEL_H

//#include <OsiSolverInterface.hpp>
#include <OsiClpSolverInterface.hpp>
#include <OsiConicSolverInterface.hpp>
#include "ConicConstraints.hpp"
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

class ColaModel: virtual public OsiConicSolverInterface,
		 public OsiClpSolverInterface {
  // data members
  ConicConstraints * cc_;
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
  // PRIVATE FUNCTIONS
  // reduce conic constraint given by size and members to smaller cones, save
  // them in reduced_cc_i
  // num_var will store the new number of variables after reduction
  void reduce_cone(int size, const int * members, ConicConstraints * reduced_cc_i, int & num_var);
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
public:
  // get problem status
  ProblemStatus problem_status();
  // functions
  ColaModel();
  ColaModel(char * data_file);
  virtual ~ColaModel();
  //void read(const char * data_file);
  // following three is useful when we clone
  void setConicConstraints(ConicConstraints * cc);
  const ConicConstraints * get_conic_constraints() const;
  void setOptions(Options * opt);
  Options * options();
  void print_stats() const;
  void print_solution() const;
  void report_feasibility() const;
  // VIRTUAL FUNCTIONS
  // get conic constraints
  virtual void getConicConstraint(int index, OsiConeType & type,
				  int & numMembers,
				  int *& members) const;
  // add conic constraints
  virtual void addConicConstraint(OsiConeType type,
				  int numMembers,
				  const int * members);
  virtual void removeConicConstraint(int index);
  //ProblemStatus solve();
  virtual int getNumCones() const;
  virtual int readMps(const char * filename, const char * extension="mps");
  virtual void initialSolve();
  virtual void resolve();
  virtual OsiConicSolverInterface * clone (bool copyData=true) const;
  // END OF VIRTUAL FUNCTIONS
  // solves problem using ben-tal nemirovski approximation
  // v is approximation parameter
  ProblemStatus solve_with_bn(bool resolve, int v);
};

#endif
