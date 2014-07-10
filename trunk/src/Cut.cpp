#include "Cut.hpp"
#include <algorithm>

Cut::Cut(const double * coef, const int * index, double rhs, int size, int cgc)
  : cut_generating_cone_(cgc), size_(size) {
  coef_ = new double[size];
  std::copy(coef, coef+size, coef_);
  index_ = new int[size];
  std::copy(index, index+size, index_);
  rhs_ = rhs;
}

// copy constructor
Cut::Cut(const Cut & other)
  : size_(other.size()),
    cut_generating_cone_(other.cut_generating_cone()) {
  // copy index
  index_ = new int[size_];
  std::copy(other.index(), other.index()+size_, index_);
  // copy coef
  coef_ = new double[size_];
  std::copy(other.coef(), other.coef()+size_, coef_);
  // copy rhs
  rhs_ = other.rhs();
}

// copy assignment operator
Cut & Cut::operator=(const Cut & rhs) {
  Cut * new_cut = new Cut(rhs);
  return *new_cut;
}


Cut::~Cut() {
  if (coef_)
    delete[] coef_;
  if (index_)
    delete[] index_;
}

int Cut::size() const {
  return size_;
}

const double * Cut::coef() const {
  return coef_;
}

double Cut::rhs() const {
  return rhs_;
}

const int * Cut::index() const {
  return index_;
}

int Cut::cut_generating_cone() const {
  return cut_generating_cone_;
}
