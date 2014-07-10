#ifndef SEPARATE_H
#define SEPARATE_H

#include "ConicConstraints.hpp"
#include "Options.hpp"
#include "Cut.hpp"

class Separate {
  std::vector<Cut*> * cuts_;
  // true if point is feasible for all cones.
  bool feasible_;
  // dimension of problem
  const int size_;
  // pointer to options
  const Options * options_;
  // return true if point is feasible for cone, else return false
  // if it is false it adds the created cut to the cuts_
  bool generate_cut(const int * cone_members, int cone_size,
		    ConeType cone_type, const double * point, int cgc);
  // No default constructor
public:
  // chec feasibility of the point, update feasible_, coef_, rhs_
  Separate(const ConicConstraints * cc, const double * point, int num_cols, const Options * options_);
  ~Separate();
  // returns true if point is feasible for all cones.
  bool is_feasible() const;
  // returns coefficient array
  const std::vector<Cut*> * cuts() const;
};

#endif
