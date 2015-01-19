#ifndef CONE_H
#define CONE_H

#include <CoinPackedVector.hpp>

typedef enum {
  LORENTZ=0,
  SCALED
} ConeType;

#define EPS 1e-5

// Abstract base class for a conic constraint
class Cone {
  ConeType type_;
public:
  Cone();
  Cone(ConeType type);
  virtual Cone * clone() const=0;
  virtual ~Cone();
  virtual ConeType type() const = 0;
  // returns 0 if point is epsilon feasible, nonzero otherwise
  //virtual int separate(int size, double const * point, int * & coef_ind,
  //		       double * & coef_val, double & rhs) const;
  virtual int separate(int size, double const * point,
		       CoinPackedVector * & cut,
		       double & rhs) const = 0;
  // size of cone, for Lorentz cones number of variables in the cone,
  // for Scaled cones number of rows of A plus 1
  virtual int size() const = 0;
  // return the feasibility of point
  // for Lorentz cones; x1-|x_2:n| or 2x1x2-|x_3:n|^2
  // for Scaled cones dx-h-|Ax-b|
  virtual double feasibility(double const * point) const = 0;
  virtual double feasibility(int size, CoinPackedVector const & point) const = 0;
};

#endif

