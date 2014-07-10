#include "ConicConstraints.hpp"
#include <algorithm>
#include <iostream>
#include <iomanip>

ConicConstraints::ConicConstraints() {
  num_cones_ = 0;
  num_members_.empty();
  members_.empty();
  type_.empty();
}

ConicConstraints * ConicConstraints::clone() const {
  ConicConstraints * new_cc = new ConicConstraints();
  for (int i=0; i<num_cones(); ++i) {
    new_cc->add_cone(cone_size(i), cone_members(i), type(i));
  }
  return new_cc;
}

ConicConstraints::~ConicConstraints() {
  for (int i=0; i<num_cones_; ++i) {
    delete[] members_[i];
  }
  num_members_.empty();
  members_.empty();
  type_.empty();
}

void ConicConstraints::add_cone(int num_members, const int * members, ConeType type) {
  if (num_members==0)
    throw EmptyCone();
  num_cones_ = num_cones_+1;
  num_members_.push_back(num_members);
  int * m = new int[num_members];
  std::copy(members, members+num_members, m);
  members_.push_back(m);
  type_.push_back(type);
}

void ConicConstraints::dump_cones() const {
  std::cout << std::setw(10) << "Cone index"
	    << std::setw(10) << "Type"
            << std::setw(10) << "Num memb."
            << std::setw(10) << "Members"
            << std::endl;
  for (int i=0; i<num_cones_; ++i) {
    std::cout << std::setw(10) << i
	      << std::setw(10) << type_[i]
              << std::setw(10) << num_members_[i]
              << std::setw(10) << members_[i][0]
              << std::endl;
    for (int j=1; j<num_members_[i]; ++j) {
      std::cout << std::setw(30) << members_[i][j] << std::endl;
    }
  }
}

void ConicConstraints::dump_cones_brief() const {
  std::cout << std::setw(10) << "Cone index"
	    << std::setw(10) << "Type"
            << std::setw(10) << "Num memb."
            << std::endl;
  for (int i=0; i<num_cones_; ++i) {
    std::cout << std::setw(10) << i
	      << std::setw(10) << type_[i]
              << std::setw(10) << num_members_[i]
              << std::endl;
  }
}

int ConicConstraints::num_cones() const {
  return num_cones_;
}

int ConicConstraints::cone_size(int i) const {
  return num_members_[i];
}

const int * ConicConstraints::cone_members(int i) const {
  return members_[i];
}

ConeType ConicConstraints::type(int i) const {
  return type_[i];
}
