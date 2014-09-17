#ifndef CONIC_CONSTRAINTS_H
#define CONIC_CONSTRAINTS_H

#include <exception>
#include <vector>

typedef enum {
  QUAD=0,
  RQUAD,
} ConeType;

class EmptyCone: std::exception {
};

class ConeDoesNotExist: std::exception {
};

class ConicConstraints {
  int num_cones_;
  std::vector<int> num_members_;
  std::vector<int *> members_;
  std::vector<ConeType> type_;
public:
  // default constructor, creates emty constraint set
  ConicConstraints();
  // copy constructor
  // copy assignment operator
  // returns a clone of this
  ConicConstraints * clone() const;
  ~ConicConstraints();
  void add_cone(int num_members, const int * members, ConeType type);
  // mainly used for debugging purposes
  void dump_cones() const;
  void dump_cones_brief() const;
  int num_cones() const;
  int cone_size(int i) const;
  const int  * cone_members(int i) const;
  ConeType type(int i) const;
  void remove_cone(int i);
  // reset data
  void reset();
};

#endif

