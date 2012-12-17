/**
 *  fld.cpp
 *  express
 *
 *  Created by Adam Roberts on 3/20/11.
 *  Copyright 2011 Adam Roberts. All rights reserved.
 */

#include "fld.h"
#include "main.h"
#include <numeric>
#include <boost/assign.hpp>
#include <iostream>
#include <fstream>
#include <boost/math/distributions/normal.hpp>

using namespace std;

const vector<double> KERNEL = boost::assign::list_of(-2.7725887222397811)
                                                    (-1.3862943611198906)
                                                    (-0.98082925301172619)
                                                    (-1.3862943611198906)
                                                    (-2.7725887222397811);

FLD::FLD(string param_file_name) {
  ifstream infile (param_file_name.c_str());
  size_t BUFF_SIZE = 99999;
  char line_buff[BUFF_SIZE];
  
  if (!infile.is_open()) {
    cerr << "ERROR: Unable to open paramater file '" << param_file_name
         << "'.\n";
    exit(1);
  }
  
  infile.getline (line_buff, BUFF_SIZE, '\n');
  infile.getline (line_buff, BUFF_SIZE, '\n');
  char *p = strtok(line_buff, "\t");
  size_t i = 0;
  
  _tot_mass = 0;
  _sum = 0;
  do {
    double val = strtod(p,NULL);
    _hist.push_back(log(val));
    _tot_mass += val;
    _sum += i*val;
    i++;
    p = strtok(NULL, "\t");
  } while (p);
  
  _tot_mass = log(_tot_mass);
  _sum = log(_sum);
  _min = max_val();;
}

FLD::FLD(double alpha, size_t max_val, size_t mean, size_t std_dev)
    : _hist(max_val+1),
      _tot_mass(LOG_0),
      _sum(LOG_0),
      _min(max_val) {
  assert(KERNEL.size() % 2 == 1);
  boost::math::normal norm(mean, std_dev);

  double tot = log(alpha*max_val);

  // Set to prior distribution
  for (size_t i = 0; i <= max_val; ++i) {
    double norm_mass = boost::math::cdf(norm,i+0.5) -
                       boost::math::cdf(norm,i-0.5);
    // Use float here so that output to the params file will be all non-zero.
    double mass = LOG_EPSILON;
    if (norm_mass != 0) {
      mass = tot + log(norm_mass);
    }
    _hist[i] = mass;
    _sum = log_add(_sum, log((double)i)+mass);
    _tot_mass = log_add(_tot_mass, mass);
  }
}

size_t FLD::max_val() const {
  return _hist.size()-1;
}

size_t FLD::min_val() const {
  if (_min == max_val()) {
      return 1;
  }
  return _min;
}

void FLD::add_val(size_t len, double mass) {
  assert(!isnan(mass));

  if (len > max_val()) {
      len = max_val();
  }
  if (len < _min) {
    _min = len;
  }

  size_t offset = len - KERNEL.size()/2;

  for (size_t i = 0; i < KERNEL.size(); i++) {
    if (offset > 0 && offset <= max_val()) {
      double k_mass = mass + KERNEL[i];
      _hist[offset] = log_add(_hist[offset], k_mass);
      _sum = log_add(_sum, log((double)offset)+k_mass);
      _tot_mass = log_add(_tot_mass, k_mass);
    }
    offset++;
  }
}

double FLD::pmf(size_t len) const {
  if (len > max_val()) {
    len = max_val();
  }
  return _hist[len]-_tot_mass;
}

vector<double> FLD::cmf() const {
  double cum = LOG_0;
  vector<double> cdf(_hist.size());
  for (size_t i = 0; i < _hist.size(); ++i) {
    cum = log_add(cum, _hist[i]);
    cdf[i] = cum - _tot_mass;
  }
  assert(approx_eq(cum, _tot_mass));

  return cdf;
}

double FLD::tot_mass() const {
    return _tot_mass;
}

double FLD::mean() const {
    return _sum - tot_mass();
}

string FLD::to_string() const {
  string s = "";
  char buffer[50];
  for(size_t i = 0; i <= max_val(); i++) {
    sprintf(buffer, "%e\t",sexp(pmf(i)));
    s += buffer;
  }
  s.erase(s.length()-1,1);
  return s;
}

void FLD::append_output(ofstream& outfile) const {
  outfile << ">Fragment Length Distribution (0-" << max_val() << ")\n";
  outfile << to_string() << endl;
}
