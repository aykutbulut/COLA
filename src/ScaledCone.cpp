#include "ScaledCone.hpp"
#include "ColaModel.hpp"

#include <numeric>
#include <cmath>

ScaledCone::ScaledCone(CoinPackedMatrix const * const A,
		       CoinPackedVector const * const b,
		       CoinPackedVector const * const d, double h)
  : Cone(SCALED) {
  A_ = new CoinPackedMatrix(*A);
  b_ = new CoinPackedVector(*b);
  d_ = new CoinPackedVector(*d);
  h_ = h;
  dense_b_ = 0;
  // compute dense_b
  compute_dense_b();
}

ScaledCone::ScaledCone(ScaledCone const & other)
  : Cone(SCALED) {
  A_ = new CoinPackedMatrix(*(other.matrixA()));
  b_ = new CoinPackedVector(*(other.vectorb()));
  d_ = new CoinPackedVector(*(other.vectord()));
  h_ = other.h();
  dense_b_ = 0;
  compute_dense_b();
}

ScaledCone & ScaledCone::operator=(ScaledCone const & rhs) {
  if (this!=&rhs) {
    Cone::operator=(rhs);
    A_ = new CoinPackedMatrix(*(rhs.matrixA()));
    b_ = new CoinPackedVector(*(rhs.vectorb()));
    d_ = new CoinPackedVector(*(rhs.vectord()));
    h_ = rhs.h();
    dense_b_ = 0;
    compute_dense_b();
  }
  return *this;
}

Cone * ScaledCone::clone() const {
  Cone * c = new ScaledCone(*this);
  return c;
}


ScaledCone::~ScaledCone() {
  delete A_;
  delete b_;
  delete d_;
  if (dense_b_)
    delete[] dense_b_;
}

CoinPackedMatrix const * ScaledCone::matrixA() const {
  return A_;
}

CoinPackedVector const * ScaledCone::vectorb() const {
  return b_;
}

CoinPackedVector const * ScaledCone::vectord() const {
  return d_;
}

double ScaledCone::h() const {
  return h_;
}

// returns 0 if point is epsilon feasible, nonzero otherwise
int ScaledCone::separate(int size, double const * point,
			 CoinPackedVector * & cut,
			 double & rhs) const {
  std::cerr << "Not implemented yet." << std::endl;
  throw std::exception();
  double feas = 0.0;
  if (feas>-options()->get_dbl_option(TOL))
    return 0;
  int n = A_->getNumCols();
  if (size!=n) {
    std::cerr << "Point size should mmatch column size of matrix A."
	      << std::endl;
    throw std::exception();
  }
  double * sol = new double[size];
  simple_separation(point, sol);
}

// size of cone, for Scaled cones number of rows of A plus 1
int ScaledCone::size() const {
  int size;
  size = A_->getNumRows()+1;
  return size;
}

// return the feasibility of point
// for Scaled cones dx-h-|Ax-b|
double ScaledCone::feasibility(double const * point) const {
  // number of rows
  int m = A_->getNumRows();
  double feas;
  double dx;
  double * Ax;
  // Ax-b
  double * Ax_b;
  // |Ax-b|
  double norm_Ax_b;
  // compute dx
  int num_elem = d_->getNumElements();
  int const * d_ind = d_->getIndices();
  double const * d_val = d_->getElements();
  dx = 0.0;
  for (int i=0; i<num_elem; ++i) {
    dx += point[d_ind[i]]*d_val[i];
  }
  // end of computing dx
  // compute Ax
  Ax = new double[m];
  A_->times(point, Ax);
  // end of computing Ax
  // compute Ax-b
  Ax_b = new double[m];
  for (int i=0; i<m; ++i) {
    Ax_b[i] = Ax[i]-dense_b_[i];
  }
  // end of computing Ax-b
  norm_Ax_b = std::inner_product(Ax_b, Ax_b+m, Ax_b, 0.0);
  norm_Ax_b = sqrt(norm_Ax_b);
  delete[] Ax;
  delete[] Ax_b;
  feas = dx-h_-norm_Ax_b;
  return feas;
}

double ScaledCone::feasibility(int size, CoinPackedVector const & point) const {
  int n = A_->getNumCols();
  double * p = new double[n]();
  int num_elem = point.getNumElements();
  int const * ind = point.getIndices();
  double const * val = point.getElements();
  for (int i=0; i<num_elem; ++i) {
    p[ind[i]] = val[i];
  }
  double feas = feasibility(p);
  delete[] p;
  return feas;
}

// go along d until you hit boundry. sol is on the boundry
// steplength alpha is solution of the following quadratic
// alpha^2(uu-(dd)^2)+4alpha(uw+hdd)+(ww-hh)=0
//
void ScaledCone::simple_separation(double const * point, double * sol) const {
  int n = A_->getNumCols();
  int m = A_->getNumRows();
  double * Ad;
  double * Ax;
  double * Ax_b;
  double dd;
  // compute Ad
  Ad = new double[m];
  A_->times(*d_, Ad);
  double * u = Ad;
  // end of computing Ad
  // compute Ax
  Ax = new double[m];
  A_->times(point, Ax);
  // end of computing Ax

  double * w = Ax_b;
  delete[] Ad;
  delete[] Ax;
  delete[] Ax_b;

}

void ScaledCone::compute_dense_b() {
  int n = A_->getNumCols();
  dense_b_ = new double[n]();
  int const * ind = b_->getIndices();
  double const * val = b_->getElements();
  int num_elem = b_->getNumElements();
  for (int i=0; i<num_elem; ++i) {
    dense_b_[ind[i]] =  val[i];
  }
}

// initial linear relaxation of conic constraints
// add dx-h>=0 for SCALED cones.
void ScaledCone::relax (ColaModel & model) const {
  model.addRow(*d_, h_, model.getInfinity());
}

// reduces conic constraint to a set of conic constraints of smaller size.
// used for bet-tal nemirovski method
std::vector<Cone*> ScaledCone::reduce() const {
  std::cerr << "Not implemented yet." << std::endl;
  throw std::exception();
}
