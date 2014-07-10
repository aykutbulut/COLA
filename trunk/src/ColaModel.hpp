#ifndef COLA_MODEL_H
#define COLA_MODEL_H

//#include <OsiSolverInterface.hpp>
#include <OsiClpSolverInterface.hpp>
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

class ColaModel: virtual public OsiClpSolverInterface {
  // data members
  ConicConstraints * cc_;
  Options * options_;
  ProblemStatus soco_status_;
  // number of cuts generated for each cone
  int * num_cuts_;
  // total number of cuts
  int total_num_cuts_;
  // number of supports generated for each cone
  int * num_supports_;
  // total number of supports
  int total_num_supports_;
  // get problem status
  ProblemStatus problem_status();
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
public:
  // functions
  ColaModel();
  ColaModel(char * data_file);
  ColaModel * clone (bool copyData=true) const;
  ~ColaModel();
  // this is for reading conic mps
  void read(const char * data_file);
  void print_stats() const;
  void print_solution() const;
  void report_feasibility() const;
  ProblemStatus solve();
  ConicConstraints * get_conic_constraints();
  Options * options();
  // set solver
  void setSolver(OsiSolverInterface * solver);
  // set conic constraints
  void setCC(ConicConstraints * cc);
  // set options
  void setOptions(Options * opt);
  //void setColBounds (const double * lb, const double * ub);
  // override some virtual functions
  virtual void initialSolve();
  virtual void resolve();
  virtual int readMps(const char * filename, const char * extension="mps");
};

#endif
