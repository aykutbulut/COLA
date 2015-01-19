#ifndef LORENTZ_CONE_H
#define LORENTZ_CONE_H

#include "Cone.hpp"

typedef enum {
  QUAD=0,
  RQUAD,
} LorentzConeType;


class LorentzCone: virtual public Cone {
  LorentzConeType lc_type_;
  int size_;
  int * members_;
  void simple_separation(double const * p, double * coef) const;
  void closest_point_separation(double const * p, double * coef) const;
  void find_closest_point(double const * y, double * sol) const;
public:
  LorentzCone(LorentzConeType type, int size, int const * members);
  // copy constructor
  LorentzCone(LorentzCone const & other);
  // copy assignment operator
  LorentzCone & operator=(LorentzCone const & rhs);
  ~LorentzCone();
  int const * members() const;
  LorentzConeType lorentz_cone_type() const;
  // VIRTUAL FUNCTIONS
  // return pointer to a clone of this
  virtual Cone * clone() const;
  virtual ConeType type() const;
  // returns 0 if point is not epsilon feasible, nonzero otherwise
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
};

#endif
