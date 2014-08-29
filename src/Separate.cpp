#include "Separate.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <numeric>
#include <algorithm>

// todo(aykut) try passing cuts_ pointer to functions, which may be more
// efficient.
Separate::Separate(const ConicConstraints * cc, const double * point,
		   int num_cols, const Options * options)
  : feasible_(true), size_(num_cols), options_(options) {
  cuts_ = new std::vector<Cut*>();
  cuts_->empty();
  const int num_cones = cc->num_cones();
  bool feas;
  for (int i=0; i<num_cones; ++i) {
    int cone_size = cc->cone_size(i);
    const int * cone_members = cc->cone_members(i);
    ConeType cone_type = cc->type(i);
    // create cut and add it to cut vector if possible
    feas = generate_cut(cone_members, cone_size, cone_type, point, i);
    if (feas) {
      continue;
    }
    feasible_ = false;
    if (options_->get_int_option(LOG_LEVEL)>0) {
      std::cout << "Cut generated using cone " << i
		<< "." << std::endl;
    }
  }
}

Separate::~Separate() {
  // free cuts_
  std::vector<Cut*>::iterator it;
  for (it=cuts_->begin(); it!=cuts_->end(); ++it) {
    delete *it;
  }
  cuts_->empty();
  delete cuts_;
}

const std::vector<Cut*> * Separate::cuts() const {
  return cuts_;
}

bool Separate::is_feasible() const {
  return feasible_;
}

// return true if point is feasible for cone, else return false
// if it is false it adds the created cut to the cuts_
bool Separate::generate_cut(const int * cone_members,
			    int cone_size,
			    ConeType cone_type,
			    const double * point,
			    int cone_index) {
  // compute coef, [2x1, -2x2, -2x3, ... -2xn]
  double * coef = new double[cone_size];
  double lhs;
  double * par_point = new double[cone_size];
  for(int i=0; i<cone_size; ++i) {
    par_point[i] = point[cone_members[i]];
  }
  double * p = par_point;
  // we want some sort of scaling of vector p with alpha.
  // note that multiplying it with alpha will not change the support we
  // will generate theoretically, but it will benefit us numerically
  // for getting cuts better in terms of numerics.
  // attempt 1, scale it to unit l_2 ball
  double sum2n = std::inner_product(p+1, p+cone_size, p+1, 0.0);
  double sqrt_sum = sqrt(sum2n+p[0]*p[0]);
  double alpha = 1.0 / sqrt_sum;
  // for (int i=0; i<cone_size; ++i) {
  //   p[i] = alpha*p[i];
  // }
  // recalculate sums
  sum2n = std::inner_product(p+1, p+cone_size, p+1, 0.0);
  if (cone_type==QUAD) {
    // set coef, lhs
    double x1 = sqrt(sum2n);
    // cone is in canonical form
    for (int i=1; i<cone_size; ++i) {
      coef[i] = 2.0*p[i];
    }
    coef[0] = -2.0*x1;
    lhs = p[0] - x1;
  }
  else {
    //  at the end, set coef and lhs
    // map point from RQUAD space to QUAD space, find the projection on QUAD,
    // project this point to RQUAD and generate cut
    double x1 = 0.0;
    double x2 = 0.0;
    // cone is a rotated cone
    // from point we move along [2point_2 2point_1 0 ... 0] until we hit
    // boundary. Then from this point in boundry we generate coef.
    // first compute u, step length
    double p1 = p[0];
    double p2 = p[1];
    x2 = (p2-p1)/2.0;
    double ssum = std::inner_product(p+2, p+cone_size, p+2, 0.0);
    double sqrt_sum = sqrt(ssum+p1*p1+p2*p2);
    double alpha = 1.0 / sqrt_sum;
    // for (int i=0; i<cone_size; ++i) {
    //   p[i] = alpha*p[i];
    // }
    p1 = p[0];
    p2 = p[1];
    // recalculate sums
    ssum = std::inner_product(p+2, p+cone_size, p+2, 0.0);
    x1 = (sqrt((-p1+p2)*(-p1+p2)+2.0*ssum) - (-p1+p2)) / 2.0;
    x2 = (sqrt((-p1+p2)*(-p1+p2)+2.0*ssum) + (-p1+p2)) / 2.0;
    // improve x1 by adjusting it a little, 2x1x2-ssum=0
    // x1 = ssum/(2.0*x2);
    // measure lhs in RQUAD space
    lhs = 2.0*p1*p2-ssum;
    // generate cut from xbar
    coef[0] = -2.0*x2;
    coef[1] = -2.0*x1;
    for (int i=2; i<cone_size; ++i) {
      coef[i] = 2.0*p[i];
    }
  }
  if (options_->get_int_option(LOG_LEVEL) > 0) {
    if (cone_type==QUAD)
      std::cout << "Cone infeasibility x_1 - sqrt(x_{2:n}^T x_{2:n}) is " << lhs << std::endl;
    else
      std::cout << "Cone infeasibility 2x_1x_2 - x_{3:n}^T x_{3:n} is " << lhs << std::endl;
  }
  // constraint is feasible if lhs is greater than -tol
  if (lhs > -options_->get_dbl_option(TOL)) {
    // point is feasible
    return true;
  }
  // if sqrt_sum is close to 0, we assume we are feasible
  if (sqrt_sum < options_->get_dbl_option(TOL)) {
    // point is feasible
    return true;
  }

  // point is not feasible, add cut to cuts_ and return false
  // index is cone_members
  // rhs is allways 0.0
  // check if we actually cut the point
  double term1 = std::inner_product(coef, coef+cone_size, p, 0.0);
  if (term1< -options_->get_dbl_option(TOL)) {
    throw "Generated plane does not cut point.";
  }
  Cut * cut = new Cut(coef, cone_members, 0.0, cone_size, cone_index);
  cuts_->push_back(cut);
  delete[] coef;
  delete[] par_point;
  return false;

}
