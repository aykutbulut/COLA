#ifndef SCALED_CONE_H
#define SCALED_CONE_H

#include "Cone.hpp"
#include "LorentzCone.hpp"
#include <CoinPackedMatrix.hpp>

// scaled cone is in |Ax-b| <= dx-h form
class ScaledCone: virtual public Cone {
  // conic constraint data
  CoinPackedMatrix * A_;
  CoinPackedVector * b_;
  CoinPackedVector * d_;
  double h_;
  // // number of columns of matrix A
  // int n_;
  // // number of rows of matrix A
  // int m_;
  // data that does not depend on x_bar, we can compute this on constructor
  double * dense_b_;
  double * dense_d_;
  // scalar d^Td
  double dd_;
  // vector Ad
  double * Ad_;
  // scalar (Ad)^T Ad
  double AdAd_;
  // matrix A^TA - dd^T, needed for computing gradient
  double * AA_dd_;
  // vector A^Tb, needed for computing gradient
  double * Ab_;
  // null space of A
  double * H_;
  // compute v, v is the smallest angle column of H, ie, max d^T h_i
  double * v_;
  // private functions
  void compute_dense_b();
  void compute_dense_d();
  void compute_dd();
  void compute_Ad();
  void compute_AdAd();
  void compute_AA_dd();
  void compute_Ab();
  void compute_H();
  void compute_Ax_b(double const * point, double * Ax_b) const;
  double compute_dx(double const * point) const;
  double compute_Ax_b_Ad(double const * point, double const * Ax_b) const;
  double compute_alpha(double const * Ax_b, double dx, double Ax_b_Ad) const;
  double compute_null_alpha(double const * Ax_b, double dx) const;
  // compute normal of tangent plane to cone at x_hat
  // stores the output at grad_x
  void compute_gradient(double const * x_hat, int & size, int *& ind_grad_x,
                        double *& val_grad_x) const;
  double compute_rhs(double const * point, int size, int const * ind,
                     double const * val) const;
  // computes step length alpha, x_hat = point + alpha d
  double compute_alpha();
  void simple_separation(double const * point, double const * Ax_b,
                         double dx, double Ax_b_Ad, double * x_hat) const;
  void null_separation(double const * point, double const * Ax_b,
                         double dx, double * x_hat) const;
  // returns a matrix with right hand size and a lorentz cone that
  // represents the scaled cone together.
  void canonical_form(CoinPackedMatrix *& mat, double & rhs,
                      LorentzCone *& c) const;
  double pos_root_of_quad_formula(double a, double b, double c) const;
  void roots(double a, double b, double c,
      double & root1, double & root2) const;
public:
  ScaledCone(CoinPackedMatrix const * A,
	     CoinPackedVector const * b,
	     CoinPackedVector const * d, double h);
  // copy constructor
  ScaledCone(ScaledCone const & other);
  // copy assignment operator
  ScaledCone & operator=(ScaledCone const & rhs);
  virtual ~ScaledCone();
  CoinPackedMatrix const * matrixA() const;
  CoinPackedVector const * vectorb() const;
  CoinPackedVector const * vectord() const;
  double h() const;
  // VIRTUAL FUNCTIONS
  // return pointer to a clone of this
  virtual Cone * clone() const;
  // returns 0 if point is epsilon feasible, nonzero otherwise
  //virtual int separate(int size, double const * point, int * & coef_ind,
  //		       double * & coef_val, double & rhs) const;
  virtual int separate(int size, double const * point,
		       CoinPackedVector * & cut,
		       double & rhs) const;
  // size of cone, for Lorentz cones number of variables in the cone,
  // for Scaled cones number of rows of A plus 1
  virtual int size() const;
  // return the feasibility of point
  // for Lorentz cones; x1-|x_2:n| or 2x1x2-|x_3:n|^2
  // for Scaled cones dx-h-|Ax-b|
  virtual double feasibility(double const * point) const;
  virtual double feasibility(int size, CoinPackedVector const & point) const;
  // initial linear relaxation of conic constraints
  // add dx-h>=0 for SCALED cones.
  virtual void relax (ColaModel & model) const;
  // reduces conic constraint to a set of conic constraints of smaller size.
  // used for bet-tal nemirovski method
  virtual std::vector<Cone*> reduce() const;
  // returns a matrix with right hand size and a lorentz cone that
  // represents the scaled cone together.

  // approximate the cone around given point.
  // If given point is in interior do nothing.
  // if it is on the boundry add support
  // We do not expect it to be infeasible for now. This may change in future
  // in case of infeasibility we will just call separate routine.
  virtual void approximate(double const * sol, OsiCuts * cuts);
};

#endif
