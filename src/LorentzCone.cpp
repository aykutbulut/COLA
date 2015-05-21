#include "LorentzCone.hpp"
#include "ColaModel.hpp"

#include <numeric>
#include <cmath>
#include <set>
#include <iomanip>

LorentzCone::LorentzCone(ConeType type, int size, int const * members)
  : Cone(type) {
  if (type==SCALED) {
    std::cerr << "Use ScaledCone class for scaled cones."
	      << std::endl;
    throw std::exception();
  }
  size_ = size;
  members_ = new int[size];
  std::copy(members, members+size, members_);
}

// copy constructor
LorentzCone::LorentzCone(LorentzCone const & other): Cone(other) {
  size_ = other.size();
  // copy members
  members_ = new int[size_];
  int const * other_members = other.members();
  std::copy(other_members, other_members+size_, members_);
}

// copy assignment operator
LorentzCone & LorentzCone::operator=(LorentzCone const & rhs) {
  if (this!=&rhs) {
    Cone::operator=(rhs);
    size_ = rhs.size();
    // copy members
    members_ = new int[size_];
    int const * rhs_members = rhs.members();
    std::copy(rhs_members, rhs_members+size_, members_);
  }
  return *this;
}

Cone * LorentzCone::clone() const {
  Cone * c = new LorentzCone(*this);
  return c;
}

LorentzCone::~LorentzCone() {
  delete[] members_;
}

int const * LorentzCone::members() const {
  return members_;
}

// VIRTUAL FUNCTIONS
// returns 0 if point is not epsilon feasible, nonzero otherwise
int LorentzCone::separate(int size, double const * point,
			  CoinPackedVector * & cut,
			  double & rhs) const {
  double feas = feasibility(point);
  if (feas>-options()->get_dbl_option(TOL))
    return 1;
  // coef array, [2x1, -2x2, -2x3, ... -2xn]
  double * coef = new double[size_];
  double * p = new double[size_];
  for (int i=0; i<size_; ++i) {
    p[i] = point[members_[i]];
  }
  // 2. compute point on cone and coefficient from the point
  simple_separation(p, coef);
  // closest_point_separation(p, coef);
  // check if we actually cut the point
  // todo(aykut) there is a problem here. we should check sol not p,
  // p is on the cone boundry.
  double term1 = std::inner_product(coef, coef+size_, p, 0.0);
  if (term1< -options()->get_dbl_option(TOL)) {
    std::cerr << "Generated plane does not cut point." << std::endl;
    throw std::exception();
  }
  // point is not feasible, add cut to cuts_ and return false
  // index is cone_members
  // rhs is allways 0.0
  cut = new CoinPackedVector(size_, members_, coef);
  delete[] coef;
  delete[] p;
  return 0;
}


// size of cone, for Lorentz cones number of variables in the cone
int LorentzCone::size() const {
  return size_;
}

// return the feasibility of point
// for Lorentz cones; x1-|x_2:n| or 2x1x2-|x_3:n|^2
double LorentzCone::feasibility(double const * point) const {
  double * par_point = new double[size_];
  for(int i=0; i<size_; ++i) {
    par_point[i] = point[members_[i]];
  }
  double * p = par_point;
  double term1;
  double term2;
  double feas;
  if (type()==LORENTZ) {
    term1 = p[0];
    term2 = std::inner_product(p+1, p+size_, p+1, 0.0);
    term2 = sqrt(term2);
  }
  else if (type()==RLORENTZ) {
    term1 = 2.0*p[0]*p[1];
    term2 = std::inner_product(p+2, p+size_, p+2, 0.0);
  }
  feas = term1-term2;
  delete[] par_point;
  return feas;
}

double LorentzCone::feasibility(int size, CoinPackedVector const & point) const {
  double * dense_point = new double[size]();
  int const * ind = point.getIndices();
  double const * val = point.getElements();
  int num_elem = point.getNumElements();
  for (int i=0; i<num_elem; ++i) {
    dense_point[ind[i]] = val[i];
  }
  double * par_point = new double[size_];
  for(int i=0; i<size_; ++i) {
    par_point[i] = dense_point[members_[i]];
  }
  delete[] dense_point;
  double * p = par_point;
  double term1;
  double term2;
  double feas;
  if (type()==LORENTZ) {
    term1 = p[0];
    term2 = std::inner_product(p+1, p+size_, p+1, 0.0);
    term2 = sqrt(term2);
  }
  else if (type()==RLORENTZ) {
    term1 = 2.0*p[0]*p[1];
    term2 = std::inner_product(p+2, p+size_, p+2, 0.0);
  }
  feas = term1-term2;
  delete[] par_point;
  return feas;
}

void LorentzCone::simple_separation(double const * p,
				 double * coef) const {
  double sum_rest;
  if (type()==LORENTZ) {
    sum_rest = std::inner_product(p+1, p+size_, p+1, 0.0);
    double x1 = sqrt(sum_rest);
    // cone is in canonical form
    for (int i=1; i<size_; ++i) {
      coef[i] = 2.0*p[i];
    }
    coef[0] = -2.0*x1;
  }
  else if (type()==RLORENTZ) {
    //  at the end, set coef and lhs
    // map point from RLORENTZ space to LORENTZ space, find the projection on LORENTZ,
    // project this point to RLORENTZ and generate cut
    sum_rest = std::inner_product(p+2, p+size_, p+2, 0.0);
    double x1 = 0.0;
    double x2 = 0.0;
    // cone is a rotated cone
    // from point we move along [2point_2 2point_1 0 ... 0] until we hit
    // boundary. Then from this point in boundry we generate coef.
    // first compute u, step length
    double p1 = p[0];
    double p2 = p[1];
    x2 = (p2-p1)/2.0;
    p1 = p[0];
    p2 = p[1];
    x1 = (sqrt((-p1+p2)*(-p1+p2)+2.0*sum_rest) - (-p1+p2)) / 2.0;
    x2 = (sqrt((-p1+p2)*(-p1+p2)+2.0*sum_rest) + (-p1+p2)) / 2.0;
    // generate cut from xbar
    coef[0] = -2.0*x2;
    coef[1] = -2.0*x1;
    for (int i=2; i<size_; ++i) {
      coef[i] = 2.0*p[i];
    }
  }
}

void LorentzCone::closest_point_separation(double const * p,
					double * coef) const {
  // compute point on boundry closest to given point
  double * sol = new double[size_]();
  find_closest_point(p, sol);
  if (type()==LORENTZ) {
    // cone is in canonical form
    for (int i=1; i<size_; ++i) {
      coef[i] = 2.0*sol[i];
    }
    coef[0] = -2.0*sol[0];
  }
  else if (type()==RLORENTZ) {
    coef[0] = -2.0*sol[1];
    coef[1] = -2.0*sol[0];
    for (int i=2; i<size_; ++i) {
      coef[i] = 2.0*sol[i];
    }
  }
  delete[] sol;
}

void LorentzCone::find_closest_point(double const * y,
			    double * sol) const {
  if (type()==RLORENTZ) {
    std::cerr << "Rotated cones are not implemented yet!" << std::endl;
    throw std::exception();
  }
  int dim = size_;
  double * start = new double[dim];
  start[0] = std::inner_product(y+1, y+dim, y+1, 0.0);
  start[0] = sqrt(start[0]);
  std::copy(y+1, y+dim, start+1);
  double u_start = -0.2;
  int iter_num = 0;
  int stop = 0;
  double * diff = new double[dim]();
  double max_diff = 0.0;
  // previous x iterate
  double * prev_x = new double[dim]();
  std::copy(start, start+dim, prev_x);
  // current x iterate
  double * x = new double[dim]();
  // previous u iterate
  double prev_u = u_start;
  // current u iterate
  double u = 0.0;
  // auxilary variables
  double A = 0.0;
  double B = 0.0;
  double prev_xi_2 = 0.0;
  double prev_x1_2 = 0.0;
  double prev_x1y1 = 0.0;
  double prev_xiyi = 0.0;
  double xi_2 = 0.0;
  // distance of solution to y
  double obj;
  double diff_u;
  // print iteration information
  std::cout << std::setw(4) << "Iter "
	    << std::setw(10) << "max_diff "
	    << std::setw(15) << "u "
	    << std::setw(15) << "x1 "
	    << std::setw(15) << "x1 -||x_2|| "
	    << std::setw(15) << "obj"
	    << std::endl;
  std::cout << std::string(78, '=') << std::endl;
  while(!stop) {
    // end of print
    // 1. update u
    // 1.1 compute A
    prev_x1_2 = prev_x[0]*prev_x[0];
    prev_xi_2 = std::inner_product(prev_x+1, prev_x+dim, prev_x+1, 0.0);
    prev_x1y1 = prev_x[0]*y[0];
    prev_xiyi = std::inner_product(prev_x+1, prev_x+dim, y+1, 0.0);
    A = (2.0*prev_x1_2)/(1-prev_u) + (2.0*prev_xi_2)/(1+prev_u);
    // 1.2 compute B
    B = prev_x1_2 - prev_xi_2 - (2.0*prev_x1y1)/(1-prev_u) + (2.0*prev_xiyi)/(1+prev_u);
    // 1.3 update u
    diff_u = B/A;
    // 2. update x_i
    for (int i=1; i<dim; ++i) {
      x[i] = (-prev_x[i]-prev_u*prev_x[i]+y[i]-diff_u*prev_x[i])/(1+prev_u) + prev_x[i];
    }
    // 3. update x1
    x[0] = (-prev_x[0]+prev_u*prev_x[0]+y[0]+diff_u*prev_x[0])/(1-prev_u) + prev_x[0];
    // 4. compute diff
    max_diff = 0.0;
    for (int i=0; i<dim; ++i) {
      double value = fabs(x[i]-prev_x[i]);
      if (max_diff<value)
	max_diff = value;
    }
    // stop if the iterates are close enough
    // if (max_diff<options()->get_dbl_option(TOL) or iter_num>20)
    //   stop = 1;
    // 5. update prev_x
    std::copy(x, x+dim, prev_x);
    u = diff_u + prev_u;
    prev_u = u;
    // 6. compute obj for reporting purposes
    obj = 0.0;
    for (int i=0; i<dim; ++i) {
      obj += (x[i]-y[i])*(x[i]-y[i]);
    }
    obj = sqrt(obj);
    iter_num++;
    // print iteration information
    xi_2 = std::inner_product(x+1, x+dim, x+1, 0.0);
    std::cout << std::setw(4) << iter_num << " "
	      << std::setw(10) << max_diff << " "
	      << std::setw(15) << u << " "
	      << std::setw(15) << x[0] << " "
	      << std::setw(15) << x[0] - sqrt(xi_2) << " "
	      << std::setw(15) << obj
	      << std::endl;

    if ((((x[0]-sqrt(xi_2))<options()->get_dbl_option(TOL)) and
	 -options()->get_dbl_option(TOL)<(x[0]-sqrt(xi_2))) or
	iter_num >40)
      stop = 1;
  }
  // store solution
  std::copy(x, x+dim, sol);
  // compute distance of solution to y.
  // report objective value
  std::cout << "Distance to y is " << obj << std::endl;
  delete[] diff;
  delete[] prev_x;
  delete[] x;
}


// initial linear relaxation of conic constraints
// add x_1>=0 for LORENTZ cones
// add x_1>=0, x_2>=0 for RLORENTZ cones
void LorentzCone::relax (ColaModel & model) const {
  if (type()==LORENTZ) {
    model.setColLower(members_[0], 0.0);
  }
  else if (type()==RLORENTZ) {
    model.setColLower(members_[0], 0.0);
    model.setColLower(members_[1], 0.0);
  }
}

// reduces conic constraint to a set of conic constraints of smaller size.
// used for bet-tal nemirovski method
std::vector<Cone*> LorentzCone::reduce() const {
  std::cerr << "Not implemented yet." << std::endl;
  throw std::exception();
}
