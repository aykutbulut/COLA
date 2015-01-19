#ifndef CONE_H
#define CONE_H

// CoinUtils headers
#include <CoinPackedVector.hpp>
#include <OsiSolverInterface.hpp>
// Cola headers
#include "Options.hpp"

typedef enum {
  LORENTZ=0, // lorentz cone, x_1 >= |x_2:n|
  RLORENTZ,  // rotated lorentz cone, 2x1x2>=|x_3:n|, x1>=0, x2>=0
  SCALED,    // scaled cone, dx-h >= |Ax-b|
} ConeType;

// Abstract base class for a conic constraint
class Cone {
  ConeType type_;
  Options * options_;
  // no default constructor
  Cone();
public:
  Cone(ConeType type);
  // copy constructor
  Cone(Cone const & other);
  // copy assignment operator
  Cone & operator=(Cone const & rhs);
  ConeType type() const {return type_;}
  Options const * options() const {return options_;}
  // VIRTUAL FUNCTIONS
  virtual Cone * clone() const=0;
  virtual ~Cone();
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
  // initial linear relaxation of conic constraints
  // add x_1>=0 for LORENTZ cones
  // add x_1>=0, x_2>=0 for RLORENTZ cones
  // add dx-h>=0 for SCALED cones.
  virtual void relax (OsiSolverInterface & model) const = 0;
  // reduces conic constraint to a set of conic constraints of smaller size.
  // used for bet-tal nemirovski method
  virtual std::vector<Cone*> reduce() const = 0;
};

#endif

