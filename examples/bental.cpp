#include <OsiConicSolverInterface.hpp>
#include <ColaModel.hpp>

int main(int argc, char ** argv) {
  OsiConicSolverInterface * si = new ColaModel();
  si->readMps(argv[1]);
  ColaModel * cola_ptr = dynamic_cast<ColaModel*>(si);
  cola_ptr->solve_reducing_cones(0);
  cola_ptr->report_feasibility();
  cola_ptr->print_stats();
  //cola_ptr->print_solution();
  return 0;
}
