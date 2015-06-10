#include <OsiConicSolverInterface.hpp>
#include <ColaModel.hpp>

// read problem and convert it to canonical form and
// write mps file.
int main (int argc, char ** argv) {
  OsiConicSolverInterface * solver;
  solver = new ColaModel();
  solver->readMps(argv[1]);
  // todo(aykut) implement writing problem to mps files
  ColaModel * canonical_solver = dynamic_cast<ColaModel*>(solver)->canonical_form();
  canonical_solver->writeMps("canonical");
  delete canonical_solver;
  return 0;
}
