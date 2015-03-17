#include "Options.hpp"

Options::Options() {
  set_defaults();
}

Options * Options::clone() const {
  Options * new_opt = new Options();
  // copy int options
  std::map<IntOptions, int>::const_iterator it1;
  for (it1=int_options_.begin(); it1!=int_options_.end(); ++it1) {
    new_opt->set_int_option(it1->first, it1->second);
  }
  // copy double options
  std::map<DblOptions, double>::const_iterator it2;
  for (it2=dbl_options_.begin(); it2!=dbl_options_.end(); ++it2) {
    new_opt->set_dbl_option(it2->first, it2->second);
  }
  // copy str options
  std::map<StrOptions, std::string>::const_iterator it3;
  for (it3=str_options_.begin(); it3!=str_options_.end(); ++it3) {
    new_opt->set_str_option(it3->first, it3->second);
  }
  return new_opt;
}

Options::~Options() {
  int_options_.clear();
  dbl_options_.clear();
  str_options_.clear();
}

void Options::set_int_option(IntOptions opt, int value) {
  int_options_[opt] = value;
}

void Options::set_dbl_option(DblOptions opt, double value) {
  dbl_options_[opt] = value;
}

void Options::set_str_option(StrOptions opt, const std::string value) {
  str_options_[opt] = value;
}

void Options::set_str_option(StrOptions opt, const char * value) {
  str_options_[opt] = std::string(value);
}

void Options::set_defaults() {
  // set default int options
  int_options_[LOG_LEVEL] = 0;
  int_options_[OSI_LOG_LEVEL] = 0;
  // set default float options
  dbl_options_[TOL] = 1e-5;
}

const int Options::get_int_option(IntOptions opt) const {
  return int_options_.at(opt);
}

const double Options::get_dbl_option(DblOptions opt) const {
  return dbl_options_.at(opt);
}

const std::string Options::get_str_option(StrOptions opt) const {
  return str_options_.at(opt);
}

