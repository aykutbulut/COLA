// creates cuts using general disjucntion cuts defined by Belotti et al.
// Here we will use parallel hyperplane disjunction (split cuts) to
// generate cuts.


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

// cut method
// size: size of the cone we use for cut generation
// type: type of the cone
// sol: x^0 we use
// disj_var: index of the disjunction variable
OsiConicCuts * GD1_cut(OsiConicSolverInterface const * solver);

int main(int argc, char ** argv) {
  OsiConicSolverInterface * si = new ColaModel();
  // read
  si->readMps(argv[1]);
  // solve
  si->initialSolve();
  // time to cut
  OsiConicCuts * cuts = GD1_cut(si);
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

// GD1 cut method, creates Julio's cut
//
//
//
//
OsiConicCuts * GD1_cut(OsiConicSolverInterface const * solver) {
  //si->writeMps("mir_before_cut");
  OsiConicCuts * cut = new OsiConicCuts();
  // compute Q
  // compute q
  // compute rho
  // compute tau
  // compute Q_tau
  // compute q_tau
  // compute rho_tau
  // represent quadric with conic constraints
  // construct cut object
  delete cut;
}

