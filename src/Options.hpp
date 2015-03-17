#ifndef OPTIONS_H
#define OPTIONS_H

#include <map>
#include <string>
typedef enum {
  LOG_LEVEL=0,
  OSI_LOG_LEVEL,
  //  end_of_int_options
} IntOptions;

typedef enum {
  TOL = 0,
  //  end_of_dbl_options
} DblOptions;

typedef enum {
  //  end_of_str_options
} StrOptions;

// if LOG_LEVEL options is 0, then no print.
// 1 basic printing

class Options {
  std::map<IntOptions, int> int_options_;
  std::map<DblOptions, double> dbl_options_;
  std::map<StrOptions, std::string> str_options_;
public:
  Options();
  Options * clone() const;
  ~Options();
  void set_int_option(IntOptions opt, int value);
  void set_dbl_option(DblOptions opt, double value);
  void set_str_option(StrOptions opt, const std::string value);
  void set_str_option(StrOptions opt, const char * value);
  // sets default options
  void set_defaults();
  const int get_int_option(IntOptions opt) const;
  const double get_dbl_option(DblOptions opt) const;
  const std::string get_str_option(StrOptions opt) const;
};

#endif
