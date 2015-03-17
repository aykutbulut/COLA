#ifndef SCALED_CONE_H
#define SCALED_CONE_H

#include "Cone.hpp"
#include <CoinPackedMatrix.hpp>

// scaled cone is in |Ax-b| <= dx-h form
class ScaledCone: virtual public Cone {
  CoinPackedMatrix * A_;
  CoinPackedVector * b_;
  CoinPackedVector * d_;
  double h_;
  double * dense_b_;
  void compute_dense_b();
  void simple_separation(double const * point, double * sol) const;
public:
  ScaledCone(CoinPackedMatrix const * A,
	     CoinPackedVector const * b,
	     CoinPackedVector const * d, double h);
  // copy constructor
  ScaledCone(ScaledCone const & other);
  // copy assignment operator
  ScaledCone & operator=(ScaledCone const & rhs);
  ~ScaledCone();
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
};

#endif
