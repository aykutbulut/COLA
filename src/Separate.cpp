#include "Separate.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <exception>

// todo(aykut) try passing cuts_ pointer to functions, which may be more
// efficient.
Separate::Separate(std::vector<Cone*> cones, double const * point,
		   int num_cols, Options const * options)
  : feasible_(true), options_(options) {
  cuts_.clear();
  rhs_.clear();
  generating_cone_.clear();
  bool feas;
  std::vector<Cone*>::const_iterator it;
  CoinPackedVector * cut;
  double rhs;
  int cone_index = 0;
  for (it=cones.begin(); it!=cones.end(); it++) {
    feas = (*it)->separate(num_cols, point, cut, rhs);
    if (!feas) {
      feasible_=false;
      cuts_.push_back(cut);
      rhs_.push_back(rhs);
      generating_cone_.push_back(cone_index);
      if (options_->get_int_option(LOG_LEVEL)>0) {
	std::cout << "Cut generated using cone " << cone_index
		  << "." << std::endl;
      }
    }
    cone_index++;
  }
}

std::vector<int> Separate::generating_cone() const {
  return generating_cone_;
}

Separate::~Separate() {
  // free cuts_
  std::vector<CoinPackedVector*>::iterator it;
  for (it=cuts_.begin(); it!=cuts_.end(); ++it) {
    delete *it;
  }
  cuts_.clear();
  rhs_.clear();
  generating_cone_.clear();
}

std::vector<CoinPackedVector*> Separate::cuts() const {
  return cuts_;
}

std::vector<double> Separate::rhs() const {
  return rhs_;
}

bool Separate::is_feasible() const {
  return feasible_;
}

