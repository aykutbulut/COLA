// numeric_solve method in Cola projects LP solutions to conic constraint
// boundry differently. To project LP solution we compute the closest point
// on the cone boundry. Then we compute the outer approximating hyperplane
// (cut) based on this point.
//
// Last time I checked this it did not improve the COLA performance
// significantly. Does not reduce the number of cuts needed to solve the problem.
//

#include <OsiConicSolverInterface.hpp>
#include <ColaModel.hpp>
#include <iomanip>

int main(int argc, char ** argv) {
  OsiConicSolverInterface * si = new ColaModel();
  si->readMps(argv[1]);
  ColaModel * cola_ptr = dynamic_cast<ColaModel*>(si);
  clock_t solution_start_time = clock();
  cola_ptr->solve_numeric();
  clock_t solution_duration = clock() - solution_start_time;
  ProblemStatus s = cola_ptr->problem_status();
  std::cout << "Is abandoned " << cola_ptr->isAbandoned() << std::endl;
  std::cout << "Is proven optimal " << cola_ptr->isProvenOptimal() << std::endl;
  std::cout << "Is proven primal infeasible " << cola_ptr->isProvenPrimalInfeasible() << std::endl;
  std::cout << "Is proven dual infeasible " << cola_ptr->isProvenDualInfeasible() << std::endl;
  std::cout << "Is primal objective limit reached " << cola_ptr->isPrimalObjectiveLimitReached() << std::endl;
  std::cout << "Is dual objective limit reached " << cola_ptr->isDualObjectiveLimitReached() << std::endl;
  std::cout << "Is iteration limit reached " << cola_ptr->isIterationLimitReached() << std::endl;
  ///// end of problem status debug
  if (s==OPTIMAL) {
    cola_ptr->report_feasibility();
    cola_ptr->print_stats();
    std::cout << "Obj value is " << std::setprecision (15) << cola_ptr->getObjValue() << std::endl;
    // report time
    std::cout << "Time spent on solve function call (in secs) " <<
      double(solution_duration) / double(CLOCKS_PER_SEC) << std::endl;
    // const double * sol = model->getColSolution();
    // for (int i=0; i<model->getNumCols(); ++i) {
    //   std::cout << sol[i] << std::endl;
    // }
  }















  return 0;
}
