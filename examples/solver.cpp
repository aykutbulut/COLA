// simple second order cone optimization solver
// that reads mosek type mps inputs.
//

#include <ColaModel.hpp>

int main(int argc, char ** argv) {
  ColaModel * m = new ColaModel();
  m->readMps(argv[1]);
  m->initialSolve();
  std::cout << "Optimal Solution is " << m->getObjValue()
            << std::endl;
  std::cout << "Number of iterations is " << m->getIterationCount()
            << std::endl;
  delete m;
  return 0;
}
