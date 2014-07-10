#ifndef CUT_H
#define CUT_H

// this is in sparse format.
class Cut {
  const int size_;
  double * coef_;
  double rhs_;
  // indices of the corresponding elements
  int * index_;
  // index of the cut generating cone
  const int cut_generating_cone_;
  // no default constructor
  Cut();
public:
  Cut(const double * coef, const int * index, double rhs, int size, int cgc);
  // copy constructor
  Cut(const Cut & other);
  // copy assignment operator
  Cut & operator=(const Cut & rhs);
  ~Cut();
  // get methods
  int size() const;
  const double * coef() const;
  double rhs() const;
  const int * index() const;
  int cut_generating_cone() const;
};

#endif
