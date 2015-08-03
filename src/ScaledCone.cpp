#include "ScaledCone.hpp"
#include "ColaModel.hpp"

#include <numeric>
#include <cmath>

extern "C" {
  #include <cblas.h>
  // lapack routines
  void dgesvd_(char *jobu, char *jobvt, int *m, int *n,
               double *a, int *lda, double *S, double *U,
               int *ldu, double *vt, int *ldvt, double *work,
               int *lwork, int *info);
}

#define EPS 1e-6


static void svDecompICL(int m, int n, double * A, double * sv);


ScaledCone::ScaledCone(CoinPackedMatrix const * const A,
		       CoinPackedVector const * const b,
		       CoinPackedVector const * const d, double h)
  : Cone(SCALED) {
  A_ = new CoinPackedMatrix(*A);
  //A_->reverseOrdering();
  b_ = new CoinPackedVector(*b);
  d_ = new CoinPackedVector(*d);
  h_ = h;
  dense_b_ = 0;
  dense_d_ = 0;
  Ad_ = 0;
  AA_dd_ = 0;
  Ab_ = 0;
  H_ = 0;
  v_ = 0;
  // compute dense_b
  compute_dense_b();
  compute_dense_d();
  compute_dd();
  compute_Ad();
  compute_AdAd();
  compute_AA_dd();
  compute_Ab();
  compute_H();
}

ScaledCone::ScaledCone(ScaledCone const & other)
  : Cone(SCALED) {
  A_ = new CoinPackedMatrix(*(other.matrixA()));
  //A_->reverseOrdering();
  b_ = new CoinPackedVector(*(other.vectorb()));
  d_ = new CoinPackedVector(*(other.vectord()));
  h_ = other.h();
  dense_b_ = 0;
  dense_d_ = 0;
  Ad_ = 0;
  AA_dd_ = 0;
  Ab_ = 0;
  H_ = 0;
  v_ = 0;
  compute_dense_b();
  compute_dense_d();
  compute_dd();
  compute_Ad();
  compute_AdAd();
  compute_AA_dd();
  compute_Ab();
  compute_H();

}

ScaledCone & ScaledCone::operator=(ScaledCone const & rhs) {
  if (this!=&rhs) {
    Cone::operator=(rhs);
    A_ = new CoinPackedMatrix(*(rhs.matrixA()));
    //A_->reverseOrdering();
    b_ = new CoinPackedVector(*(rhs.vectorb()));
    d_ = new CoinPackedVector(*(rhs.vectord()));
    h_ = rhs.h();
    dense_b_ = 0;
    dense_d_ = 0;
    Ad_ = 0;
    AA_dd_ = 0;
    Ab_ = 0;
    H_ = 0;
    v_ = 0;
    compute_dense_b();
    compute_dense_d();
    compute_dd();
    compute_Ad();
    compute_AdAd();
    compute_AA_dd();
    compute_Ab();
    compute_H();
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
  if (dense_b_) {
    delete[] dense_b_;
  }
  if (dense_d_) {
    delete[] dense_d_;
  }
  if (Ad_) {
    delete[] Ad_;
  }
  if (H_) {
    delete[] H_;
  }
  if (v_) {
    delete[] v_;
  }
  if (AA_dd_) {
    delete[] AA_dd_;
  }
  if (Ab_) {
    delete[] Ab_;
  }
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

// returns 1 if point is epsilon feasible, zero otherwise
int ScaledCone::separate(int size, double const * point,
			 CoinPackedVector * & cut,
			 double & rhs) const {
  double feas = feasibility(point);
  if (feas>-options()->get_dbl_option(TOL))
    return 1;
  int n = A_->getNumCols();
  if (size!=n) {
    std::cerr << "Point size should match column size of matrix A."
  	      << std::endl;
    throw std::exception();
  }
  int m = A_->getNumRows();
  // int the variable names x represents input point
  // variables that depend on x are local to separate function.
  // scalar d^Tx
  double dx = compute_dx(point);
  // vector Ax-b
  double * Ax_b = new double[m];
  compute_Ax_b(point, Ax_b);
  // scalar (Ax-b)^TAd
  double Ax_b_Ad = compute_Ax_b_Ad(point, Ax_b);
  // point on the surface of the cone
  double * x_hat = new double[n];
  // go along d until you hit the boundry.
  // computes x_hat_ which is on the surface of the cone.
  // simple_separation(point, Ax_b, dx, Ax_b_Ad, x_hat);
  // go along the smallest angle null space basis vector until
  // you hit the boundry
  null_separation(point, Ax_b, dx, x_hat);
  // compute gradient from xhat
  // gradient at point x is (A^TA-dd^T) x + (hd-A^Tb)
  int num_elem;
  int * grad_ind;
  double * grad_val;
  compute_gradient(x_hat, num_elem, grad_ind, grad_val);
  rhs = compute_rhs(x_hat, num_elem, grad_ind, grad_val);
  // negate gradient
  // for (int i=0; i<num_elem; ++i) {
  //   grad_val[i] = -grad_val[i];
  // }
  // rhs = -rhs;
  cut = new CoinPackedVector(num_elem, grad_ind, grad_val);
  delete[] grad_ind;
  delete[] grad_val;
  delete[] Ax_b;
  delete[] x_hat;
  return 0;
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

//
//
// another approach to implement this is following.
// Find alpha such that
// a alpha^2 + b alpha + c = 0 where
// a <-- (Ad)^T Ad - (d^Td)^2
// b <-- 2[(Ax-b)^T Ad - (d^Tx-h)d^Td]
// c <-- ||Ax-b||^2 - (d^Tx-h)^2
// this may seem complicated and computationally inefficient but it is not.
// we only need to compute terms that depend on x' at each iteration. The rest
// is fixed through iterations. And we have to compute the terms that depend on
// x' for checking feasibility anyway.
//
// The second approach seems to be computationally more efficient.
//
// sol is the point on the surface of the cone.
// how do we use the point on the surface of the cone? What is the gradient?
// gradient at point x is
// (A^TA-dd^T) x + (hd-A^Tb)
// we  may want to store vectors (A^TA-dd^T) and (hd-A^Tb)
void ScaledCone::simple_separation(double const * point, double const * Ax_b,
                                   double dx, double Ax_b_Ad,
                                   double * x_hat) const {
  // find x_hat_ on the cone surface
  double alpha = compute_alpha(Ax_b, dx, Ax_b_Ad);
  int n = A_->getNumCols();
  for (int i=0; i<n; ++i) {
    x_hat[i] = point[i] + alpha*dense_d_[i];
  }
}

// one way to implement this is computing the null space of matrix A.
// then we can keep the basis vectors of the null space. Then we need
// to find the basis vector that has the least angle with d, ie
// pick max d^T h_i. Then we solve the following for alpha
// alpha d^T h_i = ||Ax'-b||+h-d^Tx' where x' is the point to
// separate. then find point on the boundry as follows,
// x^hat = x' + alpha d
// Note that if max d^T h_i is negative ...
void ScaledCone::null_separation(double const * point, double const * Ax_b,
                                   double dx,
                                   double * x_hat) const {
  if (H_==0) {
    std::cerr << "H is null" << std::endl;
    throw std::exception();
  }
  // find x_hat_ on the cone surface
  double alpha = compute_null_alpha(Ax_b, dx);
  int n = A_->getNumCols();
  for (int i=0; i<n; ++i) {
    x_hat[i] = point[i] + alpha*v_[i];
  }
}

// compute null space basis of A
void ScaledCone::compute_H() {
  // for this we need a dense copy of A, unfortunately :(
  // we can free it once we computed the null space
  int n = A_->getNumCols();
  int m = A_->getNumRows();
  double * temp_A = new double[n*m]();
  int const * ind = A_->getIndices();
  double const * val = A_->getElements();
  // A_ is row major, we want tempA in col major
  for (int i=0; i<m; ++i) {
    int first = A_->getVectorFirst(i);
    int last = A_->getVectorLast(i);
    for (int j=first; j<last; ++j) {
      temp_A[i+ind[j]*m] = val[j];
    }
  }
  // Right hand side singular vectors of A
  double * sv = new double [n*n];
  svDecompICL(m, n, temp_A, sv);
  H_ = new double[(n-m)*n]();
  // Take the last n-m columns of V, lapack returns V^T
  for(int i=0; i<(n-m); ++i) {
    //cblas_dcopy(num_cols, (VT + num_rows_ + i), num_cols, (matH_+i*num_cols), 1);
    cblas_dcopy(n, (sv+m+i), n, (H_+i*n), 1);
  }
  // now H is col ordered and columns of it are basis for the null space.
  delete [] temp_A;
  delete [] sv;
  // compute v_, v is the smallest angle column of H, ie, max d^T h_i
  v_ = new double[n]();
  // H^T d
  double * Hd = new double [n-m]();
  cblas_dgemv(CblasColMajor, CblasTrans, n, n-m, 1.0, H_, n, dense_d_, 1, 0.0, Hd, 1);
  // determine max v_i
  double max_val = -1e-6;
  int max_ind = -1;
  for (int i=0; i<n-m; ++i) {
    if (Hd[i]>max_val) {
      max_val = Hd[i];
      max_ind = i;
    }
  }
  if (max_ind==-1) {
    std::cerr << "All innerproducts are negative!" << std::endl;
    throw std::exception();
  }
  // assume H is col ordered
  std::copy(H_+max_ind*n, H_+max_ind*n+n, v_);
  delete[] Hd;
}

static void svDecompICL(int m, int n, double * A, double * sv) {
  /* Do not compute U */
  char jobu = 'N';
  /* Compute all the rows of sv */
  char jobvt = 'A';
  int lda = m; //Leading dimension of A
  double s[m]; //Singular values of A
  //We are not computing U, so not need to assign memory
  double * u = NULL;
  int ldu = 1;
  //Right hand side singular vectors of A in sv
  int ldvt = n;
  // Working space size needed by lapack
  double worksize = 0.0;
  int lwork = -1; // Tells lapack to query only what is the space needed
  int info = 0;
  //  printf("fails Singular query\n");
  // Query what is the optimal size for the computation of the singular values
  dgesvd_(&jobu, &jobvt, &m, &n, A, &lda, s, u, &ldu,
          sv, &ldvt, &worksize, &lwork, &info);
  // Initilize the working space needed by lapack
  lwork = (int) worksize;
  double work[lwork];
  //   printf("fails Singular computing\n");
  // Here is where the actually computation is happening
  dgesvd_(&jobu, &jobvt, &m, &n, A, &lda, s, u, &ldu,
          sv, &ldvt, work, &lwork, &info);
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

void ScaledCone::compute_dense_d() {
  int n = A_->getNumCols();
  dense_d_ = new double[n]();
  int const * ind = d_->getIndices();
  double const * val = d_->getElements();
  int num_elem = d_->getNumElements();
  for (int i=0; i<num_elem; ++i) {
    dense_d_[ind[i]] =  val[i];
  }
}

void ScaledCone::compute_dd() {
  int n = A_->getNumCols();
  dd_ = std::inner_product(dense_d_, dense_d_+n, dense_d_, 0.0);
}

void ScaledCone::compute_Ad() {
  int m = A_->getNumRows();
  Ad_ = new double[m];
  A_->times(dense_d_, Ad_);
}

void ScaledCone::compute_AdAd() {
  int m = A_->getNumRows();
  AdAd_ = std::inner_product(Ad_, Ad_+m, Ad_, 0.0);
}

void ScaledCone::compute_AA_dd() {
  // AA_dd_ = new CoinPackedMatrix();
  // how to compute A^TA-dd^T
  // first compute AA_
  // A is given as col ordered.
  if (A_->isColOrdered()) {
    // complain about col ordered A.
    std::cerr << "conic constraint matrix A is column ordered." << std::endl;
    throw std::exception();
  }
  else {
    // A is row ordered.
  }
  // AA_ is also row ordered
  int n = A_->getNumCols();
  int m = A_->getNumRows();
  AA_dd_ = new double[n*n]();
  int const * ind = A_->getIndices();
  double const * val = A_->getElements();
  for (int i=0; i<m; ++i) {
    // compute a_ia_i^T for row i
    int first = A_->getVectorFirst(i);
    int last = A_->getVectorLast(i);
    for (int j=first; j<last; ++j) {
      // entry jj
      AA_dd_[ind[j]*n+ind[j]] += val[j]*val[j];
      for (int k=j+1; k<last; ++k) {
        // entry jk
        AA_dd_[ind[j]*n+ind[k]] += val[j]*val[k];
        // entry kj
        AA_dd_[ind[k]*n+ind[j]] += val[j]*val[k];
      }
    }
  }
  // subtract dd^T from AA^T
  ind = d_->getIndices();
  val = d_->getElements();
  int size = d_->getNumElements();
  for (int i=0; i<size; ++i) {
    // entry ii
    AA_dd_[ind[i]*n+ind[i]] -= val[i]*val[i];
    for (int j=i+1; j<size; ++j) {
      // entry ij
      AA_dd_[ind[i]*n+ind[j]] -= val[i]*val[j];
      // entry ji
      AA_dd_[ind[j]*n+ind[i]] -= val[i]*val[j];
    }
  }
}

void ScaledCone::compute_Ab() {
  int m = A_->getNumRows();
  Ab_ = new double[m];
  A_->transposeTimes(dense_b_, Ab_);
}

void ScaledCone::compute_Ax_b(double const * point, double * Ax_b) const {
  int m = A_->getNumRows();
  double * Ax = new double[m]();
  A_->times(point, Ax);
  //int n = A_->getNumCols();
  for (int i=0; i<m; ++i) {
    Ax_b[i] = Ax[i] - dense_b_[i];
  }
  delete[] Ax;
}

double ScaledCone::compute_dx(double const * point) const {
  int n = A_->getNumCols();
  double dx = std::inner_product(dense_d_, dense_d_+n, point, 0.0);
  return dx;
}

double ScaledCone::compute_Ax_b_Ad(double const * point,
                                   double const * Ax_b) const {
  int m = A_->getNumRows();
  double Ax_b_Ad = std::inner_product(Ax_b, Ax_b+m, Ad_, 0.0);
  return Ax_b_Ad;
}

// Find alpha such that
// a alpha^2 + b alpha + c = 0 where
// a <-- (Ad)^T Ad - (d^Td)^2
// b <-- 2[(Ax-b)^T Ad - (d^Tx-h)d^Td]
// c <-- ||Ax-b||^2 - (d^Tx-h)^2
double ScaledCone::compute_alpha(double const * Ax_b, double dx,
                                 double Ax_b_Ad) const {
  double quad_coef = AdAd_ - dd_*dd_;
  double lin_coef = 2.0*(Ax_b_Ad - (dx-h_)*dd_);
  int m = A_->getNumRows();
  double Ax_b_2 = std::inner_product(Ax_b, Ax_b+m, Ax_b, 0.0);
  double const_term = Ax_b_2 - (dx-h_)*(dx-h_);
  // find positive root of quadratic formula
  double alpha = pos_root_of_quad_formula(quad_coef, lin_coef, const_term);
  return alpha;
}

double ScaledCone::compute_null_alpha(double const * Ax_b, double dx) const {
  double alpha;
  double norm_Ax_b = 0.0;
  int m = A_->getNumRows();
  int n = A_->getNumCols();
  norm_Ax_b = std::inner_product(Ax_b, Ax_b+m, Ax_b, 0.0);
  norm_Ax_b = sqrt(norm_Ax_b);
  double dv = 0;
  dv = std::inner_product(dense_d_, dense_d_+n, v_, 0.0);
  alpha = (norm_Ax_b+h_-dx)/dv;
  return alpha;
}

// sol is the point on the surface of the cone.
// how do we use the point on the surface of the cone? What is the gradient?
// gradient at point x is
// (A^TA-dd^T) x + (hd-A^Tb)
// we  may want to store vectors (A^TA-dd^T) and (hd-A^Tb)
void ScaledCone::compute_gradient(double const * x_hat,
                                  int & size,
                                  int *& ind_grad_x,
                                  double *& val_grad_x) const {
  int n = A_->getNumCols();
  // (A^TA-dd^T)x
  double * AA_dd_x = new double[n];
  double * grad = new double[n];
  ind_grad_x = new int[n];
  val_grad_x = new double[n];
  cblas_dgemv(CblasColMajor, CblasTrans, n, n, 1.0, AA_dd_, n, x_hat, 1, 0.0, AA_dd_x, 1);
  for (int i=0; i<n; ++i ) {
    grad[i] = AA_dd_x[i] + h_*dense_d_[i]-Ab_[i];
  }
  size = 0;
  for (int i=0; i<n; ++i) {
    if (grad[i]<-1e-10 || grad[i]>1e-10) {
      ind_grad_x[size] = i;
      val_grad_x[size] = grad[i];
      size++;
    }
  }
  delete[] AA_dd_x;
  delete[] grad;
}

double ScaledCone::compute_rhs(double const * point, int size,
                               int const * ind,
                               double const * val) const {
  double rhs = 0.0;
  for (int i=0; i<size; ++i) {
    rhs += point[ind[i]]*val[i];
  }
  return rhs;
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

void ScaledCone::canonical_form(CoinPackedMatrix *& mat, double & rhs,
                                LorentzCone *& c) const {
  std::cerr << "Not implemented yet!" << std::endl;
  throw std::exception();
}

double ScaledCone::pos_root_of_quad_formula(double a, double b, double c) const {
  double root1, root2;
  roots(a, b, c, root1, root2);
  /**
   ** This considers the parallel only case where
   ** we need the max root
   **/
  return (std::max(root1, root2));
}

void ScaledCone::roots(double a, double b, double c,
		  double & root1, double & root2) const {
  //Check if it is a linear function
  //Returns the same value in the twoo rots when linear
  if(fabs(a) < EPS) {
    root1 = root2 = - c / b;
  }
  double discr = b*b - 4*a*c;
  //Check if roots are comples
  //Stops when complex roots
  if(discr < -EPS) {
    root1 = root2 = 1.0;
  }
  else {
    double numerator = 0.0;
    if (b > EPS)
      numerator = 0 - b - sqrt(discr);
    else
      numerator = 0 - b + sqrt(discr);
    /**
     ** If numerator is negative one of the roots
     ** will be infinity, this makes no sense in our case
     ** where the too roots are finite.
     ** TODO: Verify that this is not going to generate
     ** false results
     **/
    if (fabs(numerator) < EPS) {
      root1 = root2 = 0;
    }
    else {
      root1 = numerator / (2*a);
      root2 = (2*c) / numerator;
    }
  }
}

// approximate the cone around given point.
// If given point is in interior do nothing.
// if it is on the boundry add support
// We do not expect it to be infeasible for now. This may change in future
// in case of infeasibility we will just call separate routine.
void ScaledCone::approximate(double const * sol, OsiCuts * cuts) {
  throw "Not implemented yet!";
}
