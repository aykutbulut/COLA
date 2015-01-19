#include "Cone.hpp"

Cone::Cone() {
}

Cone::Cone(ConeType type): type_(type) {
  options_ = new Options();
}

// copy constructor
Cone::Cone(Cone const & other) {
  type_ = other.type();
  options_ = other.options()->clone();
}

// copy assignment operator
Cone & Cone::operator=(Cone const & rhs) {
  if (this != &rhs) {
    type_ = rhs.type();
    options_ = rhs.options()->clone();
  }
  return *this;
}

Cone::~Cone() {
  if (options_)
    delete options_;
}
