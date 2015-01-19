#ifndef SEPARATE_H
#define SEPARATE_H

#include "Cone.hpp"
#include "Options.hpp"
#include <vector>
#include <CoinPackedVector.hpp>

class Separate {
  std::vector<CoinPackedVector*> cuts_;
  std::vector<double> rhs_;
  std::vector<int> generating_cone_;
  bool feasible_;
  Options const * options_;
public:
  // chec feasibility of the point, update feasible_, coef_, rhs_
  Separate(std::vector<Cone*> cones, double const * point, int num_cols,
	   Options const * options);
  ~Separate();
  // returns true if point is feasible for all cones.
  bool is_feasible() const;
  // returns coefficient array
  std::vector<CoinPackedVector*> cuts() const;
  std::vector<double> rhs() const;
  std::vector<int> generating_cone() const;
};

#endif
